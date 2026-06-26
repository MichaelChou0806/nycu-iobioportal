// =====================================================================
// analyses/groupCompare.js — 「分組比較」分析模組
// ---------------------------------------------------------------------
// 一個分析模組 = 一個自給自足的檔：它擁有自己的控制項與結果區，
// 只「呼叫」core 的引擎（dataset / stats / plots），不依賴其他分析。
// 未來新增分析就照這個樣子再寫一個檔，並在 app.js 註冊一行。
//
// 對外只暴露一個物件：{ id, name, mount(container, ctx) }
// =====================================================================

import { mannWhitney, median, mean, sd, pStars } from "../core/stats.js";
import { scatterSVG, barSVG } from "../core/plots.js";
import { downloadText } from "../core/data.js";

// 可用的分組（臨床欄位 + 兩個組別值）
const GROUPINGS = {
  vital: { col: "vital_status", a: "Alive", b: "Dead", label: "Alive_vs_Dead" },
  node:  { col: "node_status",  a: "N-",    b: "N+",   label: "Nneg_vs_Npos" },
};

export const groupCompare = {
  id: "groupCompare",
  name: "Group Comparison",

  mount(container, { dataset }) {
    container.innerHTML = `
      <div class="card">
        <div class="row">
          <div class="field">
            <label>Gene (symbol or ENSG)</label>
            <input id="gc-gene" placeholder="e.g. TP53" value="TP53" autocomplete="off">
          </div>
          <div class="field">
            <label>Cancer</label>
            <select id="gc-cancer"></select>
          </div>
          <div class="field">
            <label>Grouping</label>
            <select id="gc-grouping">
              <option value="vital">Alive vs Dead (vital_status)</option>
              <option value="node">N- vs N+ (node_status)</option>
            </select>
          </div>
          <button id="gc-run">Analyze</button>
        </div>
        <div id="gc-variant" class="variant hidden"></div>
        <div id="gc-status" class="status"></div>
      </div>
      <div id="gc-result" class="card hidden">
        <div id="gc-title" class="result-title"></div>
        <div id="gc-stats" class="stats"></div>
        <div class="charts">
          <div class="chart-box"><div id="gc-chart1"></div></div>
          <div class="chart-box"><div id="gc-chart2"></div></div>
        </div>
        <div style="margin-top:18px">
          <button id="gc-export" class="ghost">Download CSV for Prism</button>
          <div class="hint">Two columns (one per group), values stacked — paste into a Prism Column table.</div>
        </div>
      </div>`;

    const $ = sel => container.querySelector(sel);

    // 填癌種下拉
    const csel = $("#gc-cancer");
    dataset.cancers.forEach(c => {
      const o = document.createElement("option");
      o.value = c.code; o.textContent = `${c.code} (tumor ${c.n_tumor})`;
      csel.appendChild(o);
    });
    if ([...csel.options].some(o => o.value === "HNSC")) csel.value = "HNSC";

    let lastResult = null;

    async function analyze() {
      const variantBox = $("#gc-variant");
      const statusEl = $("#gc-status");
      $("#gc-result").classList.add("hidden");
      statusEl.className = "status";

      const cancer = $("#gc-cancer").value;
      const G = GROUPINGS[$("#gc-grouping").value];

      // 解析基因（含重名候選處理）
      const r = dataset.resolveGene($("#gc-gene").value);
      if (r.error) { statusEl.className = "status err"; statusEl.textContent = r.error; return; }
      let rec = r.rec;
      if (r.multiple) {
        variantBox.classList.remove("hidden");
        variantBox.innerHTML = `此 symbol 對應多個基因，請選一個再按 Analyze：<br>` +
          r.multiple.map(m => `<label style="display:inline-block;margin:6px 12px 0 0">
            <input type="radio" name="gc-vpick" value="${m.gene_id}"> ${m.gene_id} · ${m.biotype}</label>`).join("");
        variantBox.querySelectorAll("input").forEach(inp =>
          inp.addEventListener("change", () => { variantBox.dataset.pick = inp.value; }));
        if (variantBox.dataset.pick) {
          rec = dataset.geneIndex.by_ensembl[variantBox.dataset.pick.split(".")[0]];
        } else { statusEl.textContent = "請先在上方選擇要分析的基因。"; return; }
      } else {
        variantBox.classList.add("hidden"); variantBox.dataset.pick = "";
      }

      statusEl.textContent = "Fetching gene data…";
      let vals;
      try { vals = await dataset.getGeneValues(rec); }
      catch (e) { statusEl.className = "status err"; statusEl.textContent = "讀取基因檔失敗（可能 CORS 未設）: " + e.message; return; }

      // join：選定癌種、只取 tumor、每病人一個、排除 redacted、排除分組值缺失
      const samples = dataset.samples, clinical = dataset.clinical;
      const seen = new Set();
      const A = [], B = [];
      for (let i = 0; i < samples.length; i++) {
        const s = samples[i];
        if (s.cancer !== cancer || s.sample_class !== "tumor") continue;
        if (seen.has(s.patient_id)) continue;
        seen.add(s.patient_id);
        const c = clinical.get(s.patient_id);
        if (!c) continue;
        if (c.redaction && c.redaction.trim() !== "") continue;
        const ev = vals[i];
        if (!isFinite(ev)) continue;     // miRNA 缺值樣本：跳過（避免污染統計）
        const gv = (c[G.col] || "").trim();
        if (gv === G.a) A.push(ev);
        else if (gv === G.b) B.push(ev);
      }

      if (A.length < 1 || B.length < 1) {
        statusEl.className = "status err";
        statusEl.textContent = `分組後其中一組沒有樣本（${G.a}=${A.length}, ${G.b}=${B.length}）。此癌種可能沒有這個臨床分組。`;
        return;
      }

      const test = mannWhitney(A, B);
      const st = {
        a: { n: A.length, med: median(A), mean: mean(A), sd: sd(A) },
        b: { n: B.length, med: median(B), mean: mean(B), sd: sd(B) },
      };
      lastResult = { symbol: rec.symbol || rec.gene_id, cancer, label: G.label, names: [G.a, G.b], values: [A, B] };

      render(rec, cancer, G, st, test, A, B);
      statusEl.textContent = "Done.";
    }

    function render(rec, cancer, G, st, test, A, B) {
      const unit = rec.assay === "mirna_rpm" ? "RPM" : "TPM";   // miRNA 標 RPM，其餘 TPM
      $("#gc-result").classList.remove("hidden");
      $("#gc-title").textContent = `${rec.symbol || rec.gene_id} expression in ${dataset.name} · ${cancer} — ${G.a} vs ${G.b}`;
      const f = x => x.toFixed(2);
      $("#gc-stats").innerHTML = `
        <div class="stat"><div class="k">${G.a} (n)</div><div class="v">${st.a.n}</div></div>
        <div class="stat"><div class="k">${G.b} (n)</div><div class="v">${st.b.n}</div></div>
        <div class="stat"><div class="k">Median ${G.a}/${G.b}</div><div class="v">${f(st.a.med)} / ${f(st.b.med)}</div></div>
        <div class="stat"><div class="k">Mean ${G.a}/${G.b}</div><div class="v">${f(st.a.mean)} / ${f(st.b.mean)}</div></div>
        <div class="stat"><div class="k">Mann–Whitney p</div><div class="v pval">${test.p < 1e-4 ? "< 0.0001" : test.p.toPrecision(3)} <span style="color:var(--muted)">${pStars(test.p)}</span></div></div>`;

      const all = A.concat(B);
      let yMax = Math.max(...all, st.a.mean + st.a.sd, st.b.mean + st.b.sd);
      if (yMax <= 0) yMax = 1; yMax *= 1.08;

      $("#gc-chart1").innerHTML = scatterSVG(
        [{ name: G.a, vals: A, med: st.a.med, color: "var(--groupA)" },
         { name: G.b, vals: B, med: st.b.med, color: "var(--groupB)" }], yMax, `${unit} (jitter + median)`);
      $("#gc-chart2").innerHTML = barSVG(
        [{ name: G.a, m: st.a.mean, sd: st.a.sd, n: st.a.n, color: "var(--groupA)" },
         { name: G.b, m: st.b.mean, sd: st.b.sd, n: st.b.n, color: "var(--groupB)" }], yMax, `${unit} (mean ± SD)`);
    }

    function exportCSV() {
      if (!lastResult) return;
      const { names, values } = lastResult;
      const maxLen = Math.max(values[0].length, values[1].length);
      let csv = names.join(",") + "\n";
      for (let i = 0; i < maxLen; i++) {
        const a = i < values[0].length ? values[0][i].toFixed(4) : "";
        const b = i < values[1].length ? values[1][i].toFixed(4) : "";
        csv += a + "," + b + "\n";
      }
      downloadText(`${lastResult.symbol}_${lastResult.cancer}_${lastResult.label}.csv`, csv);
    }

    $("#gc-run").addEventListener("click", analyze);
    $("#gc-gene").addEventListener("keydown", e => { if (e.key === "Enter") analyze(); });
    $("#gc-export").addEventListener("click", exportCSV);
  },
};
