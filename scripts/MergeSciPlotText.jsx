/**
 * MergeSciPlotText.jsx — Adobe Illustrator ExtendScript (ES3 compatible, AI CC 2018–2026 / 29.5+)
 *
 * One-click fix for fragmented scientific plot labels imported from Python (matplotlib/scanpy)
 * or R (ggplot2/cairo_pdf/ComplexHeatmap) PDFs and SVGs!
 *
 * Merges into ONE single native editable Illustrator PointText object:
 *   - Inline italic / Greek / paren splits: "Validation (" + "n" + " = 18)" -> "Validation (n = 18)"
 *                                           "ρ" + " = 0.68, " + "p" + " = 0.002"
 *   - Subscripts & superscripts:            "log" + "2" + "(Fold change)"  -> "log₂(Fold change)"
 *                                           "−log" + "10" + "(" + "p" + "adj" + ")" (0° & 90° vertical!)
 *                                           "p = 1.4 × 10" + "−9"
 *
 * Usage:
 *   1. Open your PDF or SVG in Adobe Illustrator (File -> Open; if Placed, click Embed first).
 *   2. Optionally select specific text frames (or leave unselected to process the whole document).
 *   3. Run: File -> Scripts -> Other Script... (Ctrl+F12) -> MergeSciPlotText.jsx
 */

#target illustrator

(function () {
    if (app.documents.length === 0) {
        alert("Please open a document first.");
        return;
    }

    var doc = app.activeDocument;

    function isEditableItem(item) {
        var cur = item;
        while (cur && cur.typename !== "Document") {
            try {
                if (cur.locked === true) return false;
                if (cur.hidden === true) return false;
                if (cur.visible === false) return false;
            } catch (e) {}
            try {
                cur = cur.parent;
            } catch (e) {
                break;
            }
        }
        return true;
    }

    function trimTrailingBreaks(s) {
        return s.replace(/[\r\n]+$/, "");
    }

    function isCandidateFrame(tf) {
        if (!tf || tf.typename !== "TextFrame") return false;
        if (!isEditableItem(tf)) return false;
        var txt = "";
        try { txt = tf.contents; } catch (e) { return false; }
        if (!txt || txt.length === 0) return false;
        var clean = trimTrailingBreaks(txt);
        if (clean.length === 0) return false;
        // Skip genuine multi-line paragraphs
        if (clean.indexOf("\r") !== -1 || clean.indexOf("\n") !== -1) return false;
        return true;
    }

    // Collect from selection first (if >= 2 frames selected), else scan ALL doc.textFrames
    // (doc.textFrames automatically traverses ALL sublayers, groups, and clipping masks!)
    function collectFromSelection(sel, out) {
        var i;
        for (i = 0; i < sel.length; i++) {
            var it = sel[i];
            if (it.typename === "TextFrame") {
                if (isCandidateFrame(it)) out.push(it);
            } else if (it.typename === "GroupItem") {
                collectFromSelection(it.pageItems, out);
            }
        }
    }

    var targets = [];
    if (doc.selection && doc.selection.length > 0 && doc.selection.typename !== "TextRange") {
        try {
            collectFromSelection(doc.selection, targets);
        } catch (e) {}
    }

    // Fallback to ALL editable textFrames in the entire document (handles nested sublayers & clip groups!)
    if (targets.length < 2) {
        targets = [];
        var allTF = doc.textFrames;
        for (var ti = 0; ti < allTF.length; ti++) {
            if (isCandidateFrame(allTF[ti])) {
                targets.push(allTF[ti]);
            }
        }
    }

    if (targets.length < 2) {
        var msg = "未找到足够的可编辑文字对象（当前找到 " + targets.length + " 个，总 textFrames = " + doc.textFrames.length + "）。";
        if (doc.placedItems && doc.placedItems.length > 0) {
            msg += "\n\n检测到文档中有置入的链接对象（PlacedItem），请先选中图片点顶部栏的【嵌入】，或直接使用【文件 → 打开】打开 PDF，再运行脚本。";
        } else {
            msg += "\n\n请确认图层未锁定，且文字未被转为轮廓（Create Outlines）。";
        }
        alert(msg);
        return;
    }

    // Extract geometry in local text-reading coordinates (handles 0°, ±90°, 30°, 45° and scaled groups!)
    function getFrameInfo(tf) {
        var m = tf.matrix;
        var scaleX = Math.sqrt(m.mValueA * m.mValueA + m.mValueB * m.mValueB);
        if (scaleX < 1e-6) scaleX = 1.0;
        var cosA = m.mValueA / scaleX;
        var sinA = m.mValueB / scaleX;
        var angle = Math.atan2(sinA, cosA);

        var ax = tf.anchor[0];
        var ay = tf.anchor[1];

        // u = coordinate along reading direction; v = perpendicular to baseline (upward normal)
        var u0 = ax * cosA + ay * sinA;
        var v0 = -ax * sinA + ay * cosA;

        // Visual bounds projected onto reading axis u
        var gb = tf.geometricBounds; // [left, top, right, bottom]
        var corners = [
            [gb[0], gb[1]],
            [gb[2], gb[1]],
            [gb[2], gb[3]],
            [gb[0], gb[3]]
        ];
        var uMin = 1e12, uMax = -1e12;
        for (var ci = 0; ci < 4; ci++) {
            var cu = corners[ci][0] * cosA + corners[ci][1] * sinA;
            if (cu < uMin) uMin = cu;
            if (cu > uMax) uMax = cu;
        }

        var maxSz = 0;
        var chars = tf.characters;
        for (var k = 0; k < chars.length; k++) {
            var sz = chars[k].characterAttributes.size;
            if (sz > maxSz) maxSz = sz;
        }
        if (maxSz <= 0) maxSz = 6.0;
        var effSz = maxSz * scaleX;

        // Account for leading/trailing spaces whose geometricBounds are zero-width
        var txt = trimTrailingBreaks(tf.contents);
        var leadSp = 0, trailSp = 0;
        while (leadSp < txt.length && (txt.charAt(leadSp) === " " || txt.charAt(leadSp) === "\u00a0")) leadSp++;
        while (trailSp < txt.length && (txt.charAt(txt.length - 1 - trailSp) === " " || txt.charAt(txt.length - 1 - trailSp) === "\u00a0")) trailSp++;
        var spaceW = effSz * 0.28;

        var uStart = Math.min(u0, uMin) - leadSp * spaceW;
        // If frame is pure whitespace, estimate width directly from space count
        var uEnd;
        if (leadSp === txt.length) {
            uStart = u0;
            uEnd = u0 + txt.length * spaceW;
        } else {
            uEnd = uMax + trailSp * spaceW;
            if (uStart < uMin - spaceW * (leadSp + 1.5)) {
                uStart = uMin - leadSp * spaceW;
            }
        }

        return {
            tf: tf,
            angle: angle,
            scale: scaleX,
            uStart: uStart,
            uEnd: uEnd,
            v: v0,
            size: effSz,
            rawSize: maxSz
        };
    }

    var infos = [];
    for (var i = 0; i < targets.length; i++) {
        try {
            infos.push(getFrameInfo(targets[i]));
        } catch (e) {}
    }

    function anglesMatch(a1, a2) {
        var diff = Math.abs(a1 - a2);
        while (diff > Math.PI) diff = Math.abs(diff - 2 * Math.PI);
        return diff < (2.5 * Math.PI / 180.0);
    }

    // BFS clustering of adjacent single-line fragments (including subscripts & superscripts!)
    var n = infos.length;
    var visited = [];
    for (var vi = 0; vi < n; vi++) visited[vi] = false;

    var groups = [];
    for (var i = 0; i < n; i++) {
        if (visited[i]) continue;
        var queue = [i];
        var comp = [i];
        visited[i] = true;

        while (queue.length > 0) {
            var curIdx = queue.shift();
            var A = infos[curIdx];
            for (var j = 0; j < n; j++) {
                if (visited[j]) continue;
                var B = infos[j];
                if (!anglesMatch(A.angle, B.angle)) continue;

                var refSz = Math.max(A.size, B.size);
                var dv = Math.abs(A.v - B.v);

                // Baseline tolerance:
                //   Same baseline (italic/Greek/paren switches): dv <= 0.22 * refSz
                //   Subscript / superscript (smaller fragment <= 0.88 * refSz): dv <= 0.85 * refSz
                var oneSmaller = (Math.min(A.size, B.size) <= 0.88 * refSz);
                var maxDv = oneSmaller ? (0.85 * refSz) : (0.22 * refSz);
                if (dv > maxDv) continue;

                // Gap along reading axis u
                var gap = (A.uStart <= B.uStart) ? (B.uStart - A.uEnd) : (A.uStart - B.uEnd);

                // Tight horizontal gap prevents gluing independent tick numbers (-4, -2, 0, 2, 4)!
                var maxGap = 0.55 * refSz;
                var minGap = -0.50 * refSz;
                if (gap >= minGap && gap <= maxGap) {
                    visited[j] = true;
                    queue.push(j);
                    comp.push(j);
                }
            }
        }

        if (comp.length >= 2) {
            var items = [];
            for (var ci = 0; ci < comp.length; ci++) items.push(infos[comp[ci]]);
            items.sort(function (x, y) { return x.uStart - y.uStart; });
            groups.push(items);
        }
    }

    function snapshotChars(tf, baseShiftDelta) {
        var out = [];
        var chars = tf.characters;
        for (var k = 0; k < chars.length; k++) {
            var ch = chars[k].contents;
            if (ch === "\r" || ch === "\n") continue;
            var ca = chars[k].characterAttributes;
            var fillCol = null;
            try { fillCol = ca.fillColor; } catch (e) {}
            out.push({
                ch: ch,
                font: ca.textFont,
                size: ca.size,
                baselineShift: (ca.baselineShift || 0) + baseShiftDelta,
                fillColor: fillCol
            });
        }
        return out;
    }

    var mergedGroups = 0;
    var removedFrames = 0;

    for (var gi = 0; gi < groups.length; gi++) {
        var grp = groups[gi];

        // Main baseline v = largest-size fragment's baseline
        var mainSize = 0;
        var mainV = grp[0].v;
        var mainScale = grp[0].scale;
        for (var m = 0; m < grp.length; m++) {
            if (grp[m].size > mainSize + 0.05) {
                mainSize = grp[m].size;
                mainV = grp[m].v;
                mainScale = grp[m].scale;
            }
        }

        var allChars = [];
        var fullStr = "";
        for (var m = 0; m < grp.length; m++) {
            var item = grp[m];
            // Convert document-space baseline delta back into local unscaled pt for baselineShift
            var dyShift = (item.v - mainV) / (mainScale > 1e-6 ? mainScale : 1.0);
            if (Math.abs(dyShift) < 0.15) dyShift = 0;
            var cList = snapshotChars(item.tf, dyShift);
            for (var k = 0; k < cList.length; k++) {
                allChars.push(cList[k]);
                fullStr += cList[k].ch;
            }
        }

        var host = grp[0].tf;
        // If host itself was a sub/superscript, shift host anchor onto main baseline before rebuilding
        var hostDy = grp[0].v - mainV;
        if (Math.abs(hostDy) >= 0.15) {
            var sinA = Math.sin(grp[0].angle);
            var cosA = Math.cos(grp[0].angle);
            host.translate(hostDy * sinA, -hostDy * cosA);
        }

        host.contents = fullStr;
        var hChars = host.characters;
        for (var k = 0; k < allChars.length && k < hChars.length; k++) {
            var spec = allChars[k];
            var dca = hChars[k].characterAttributes;
            try { dca.textFont = spec.font; } catch (e) {}
            try { dca.size = spec.size; } catch (e) {}
            try { dca.baselineShift = spec.baselineShift; } catch (e) {}
            if (spec.fillColor) {
                try { dca.fillColor = spec.fillColor; } catch (e) {}
            }
        }

        for (var m = 1; m < grp.length; m++) {
            try {
                grp[m].tf.remove();
                removedFrames++;
            } catch (e) {}
        }
        mergedGroups++;
    }

    alert("合并完成！共扫描 " + targets.length + " 个文字对象，将 " + (removedFrames + mergedGroups) +
          " 个碎片合并为 " + mergedGroups + " 个完整单行文本框（保留斜体与上下标）。");
})();
