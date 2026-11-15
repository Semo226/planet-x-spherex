/* Planet X and SPHEREx — epoch animation viewer (vanilla JS).
   Talks to the FastAPI backend at /api/*. */

const $ = (id) => document.getElementById(id);

let frames = [];        // original frames
let diffFrames = [];    // processed difference frames
let current = 0;
let playing = false;
let timer = null;
let viewMode = "orig";          // 'orig' | 'diff'
let diffReady = false;

function activeFrames() {
  return viewMode === "diff" ? diffFrames : frames;
}

function setStatus(msg, cls = "") {
  const el = $("status");
  el.textContent = msg;
  el.className = "status " + cls;
}

function setViewMode(mode) {
  if (mode === "diff" && !diffReady) return;
  viewMode = mode;
  $("view-orig").classList.toggle("active", mode === "orig");
  $("view-diff").classList.toggle("active", mode === "diff");
  const list = activeFrames();
  $("slider").max = Math.max(0, list.length - 1);
  showFrame(Math.min(current, list.length - 1));
}

function showFrame(i) {
  const list = activeFrames();
  if (!list.length) return;
  current = Math.max(0, Math.min(i, list.length - 1));
  const f = list[current];
  $("frame-img").src = f.png;
  $("frame-img").style.display = "block";
  $("empty").style.display = "none";
  $("m-index").textContent = `${current + 1} / ${list.length}`;
  $("m-time").textContent = f.t_min_iso || `MJD ${f.t_min?.toFixed(3) ?? "?"}`;
  $("m-obs").textContent = f.obs_id || "—";
  $("m-det").textContent = f.detector || "—";
  $("m-band").textContent = f.bandpass || "—";
  $("slider").value = current;
  $("frame-count").textContent = `${current + 1} / ${list.length}`;
}

function stop() {
  playing = false;
  $("play").textContent = "▶ Play";
  if (timer) { clearTimeout(timer); timer = null; }
}

function schedule() {
  if (!playing) return;
  const speed = parseInt($("speed").value, 10) || 1000;
  timer = setTimeout(() => {
    if (!playing) return;
    const list = activeFrames();
    if (current >= list.length - 1) { stop(); return; }
    showFrame(current + 1);
    schedule();
  }, speed);
}

async function load() {
  const ra = parseFloat($("ra").value);
  const dec = parseFloat($("dec").value);
  const fov = parseFloat($("fov").value);
  const maxFrames = parseInt($("maxframes").value, 10) || 60;
  if (isNaN(ra) || isNaN(dec)) { setStatus("Invalid coordinates.", "error"); return; }

  stop();
  frames = [];
  diffFrames = [];
  diffReady = false;
  current = 0;
  $("slider").disabled = true;
  $("play").disabled = true;
  $("speed").disabled = true;
  $("view-diff").disabled = true;
  $("view-orig").classList.add("active");
  $("view-diff").classList.remove("active");
  viewMode = "orig";
  $("frame-img").style.display = "none";
  $("empty").style.display = "block";
  $("empty").innerHTML = "Querying IRSA & downloading cutouts…<br/><small>This can take a minute.</small>";
  setStatus(`Searching SPHEREx DR1/QR3 at RA=${ra}, Dec=${dec} …`);

  try {
    const res = await fetch("/api/query", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ ra, dec, fov_size: fov, max_frames: maxFrames }),
    });
    const data = await res.json();
    if (!res.ok) {
      setStatus(`Query failed: ${data.detail || res.status}`, "error");
      $("empty").innerHTML = "Query failed.";
      return;
    }
    frames = data.frames || [];
    if (!frames.length) {
      setStatus(`No SPHEREx observations found at RA=${ra}, Dec=${dec}.`, "ok");
      $("empty").innerHTML = "No coverage at this position.<br/>Try another coordinate.";
      return;
    }
    $("m-epochs").textContent = `${data.n_epochs} epochs / ${data.n_results} MEFs`;
    $("slider").max = frames.length - 1;
    $("slider").disabled = false;
    $("play").disabled = false;
    $("speed").disabled = false;
    showFrame(0);
    const offlineNote = data.offline ? " ⚠ OFFLINE cached data (IRSA unreachable)." : "";
    setStatus(
      `Loaded ${frames.length} frames (${data.n_epochs} epochs) at RA=${data.ra}, Dec=${data.dec}.${offlineNote} Building difference stack…`,
      data.offline ? "error" : "ok"
    );
    buildDiffStack();
    lookupCatalog(data.ra, data.dec);
  } catch (e) {
    setStatus(`Request error: ${e}`, "error");
    $("empty").innerHTML = "Request failed — is the backend running?";
  }
}

function esc(s) {
  return String(s ?? "").replace(/[&<>"']/g,
    (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}

function setCatStatus(text, cls) {
  $("cat-status-text").textContent = text;
  $("cat-status").className = "cat-status " + cls;
  $("cat-spinner").hidden = cls !== "pending";
}

function renderCatalog(data) {
  const cls = data.classification === "known_object" ? "known"
    : data.classification === "discovery_candidate" ||
      data.classification === "unmatched_within_tolerance" ? "candidate"
    : "error";
  const label = {
    known_object: "Known Object",
    unmatched_within_tolerance: "Candidate (near misses)",
    discovery_candidate: "Discovery Candidate",
  }[data.classification] || "Lookup issue";
  setCatStatus(label, cls);
  $("cat-summary").textContent = data.summary || "";

  const matches = (data.matches || []).map((m) => {
    const props = [];
    if (m.g_mean_mag != null) props.push(`G=${m.g_mean_mag.toFixed(2)}`);
    if (m.w1mpro != null) props.push(`W1=${m.w1mpro.toFixed(2)}`);
    if (m.w2mpro != null) props.push(`W2=${m.w2mpro.toFixed(2)}`);
    if (m.pmra_masyr != null || m.pmdec_masyr != null)
      props.push(`pm=${(m.pmra_masyr ?? 0).toFixed(1)},${(m.pmdec_masyr ?? 0).toFixed(1)} mas/yr`);
    if (m.parallax_mas != null) props.push(`ϖ=${m.parallax_mas.toFixed(2)} mas`);
    const fl = (m.flags || []).map((f) => `<div class="cat-flag">⚠ ${esc(f)}</div>`).join("");
    return `<div class="cat-match">
      <div class="cat-name">${esc(m.catalog)} · ${esc(m.source_id)}</div>
      <div class="cat-detail">sep ${m.separation_arcsec.toFixed(2)}″ · ${props.join(" · ") || "no photometry"}</div>
      ${fl}</div>`;
  });
  $("cat-matches").innerHTML = matches.join("") ||
    '<div class="cat-detail">No sources in search radius.</div>';

  $("cat-flags").innerHTML = (data.flags || [])
    .map((f) => `<div class="cat-flag">⚠ ${esc(f)}</div>`).join("");
  $("cat-warnings").innerHTML = (data.warnings || [])
    .map((w) => `<div class="cat-warn">⚠ ${esc(w)}</div>`).join("");
  $("cat-disclaimer").textContent = data.disclaimer || "";
}

async function lookupCatalog(ra, dec) {
  setCatStatus("Querying Gaia DR3 + AllWISE…", "pending");
  $("cat-matches").innerHTML = "";
  $("cat-flags").innerHTML = "";
  $("cat-warnings").innerHTML = "";
  $("cat-disclaimer").textContent = "";
  try {
    const res = await fetch("/api/catalog", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ ra, dec, radius_arcsec: 10, epoch: 2025.5 }),
    });
    const data = await res.json();
    if (!res.ok) {
      setCatStatus("Catalog lookup failed", "error");
      $("cat-warnings").innerHTML =
        `<div class="cat-warn">⚠ ${esc(data.detail || res.status)}</div>`;
      return;
    }
    renderCatalog(data);
  } catch (e) {
    setCatStatus("Catalog lookup failed", "error");
    $("cat-warnings").innerHTML =
      `<div class="cat-warn">⚠ ${esc(e)}</div>`;
  }
}

async function buildDiffStack() {
  try {
    const res = await fetch("/api/process", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ align: true, method: "median" }),
    });
    const data = await res.json();
    if (!res.ok) {
      setStatus(`Difference imaging failed: ${data.detail || res.status}`, "error");
      return;
    }
    diffFrames = data.frames || [];
    diffReady = diffFrames.length > 0;
    $("view-diff").disabled = !diffReady;
    $("view-diff").textContent = diffReady
      ? `Processed (Difference · ${data.detector})`
      : "Processed (Difference)";
    setStatus(
      `Loaded ${frames.length} frames (${"?"} epochs). ` +
      `Difference stack ready: detector ${data.detector}, ${data.n_frames} residual frames ` +
      `(limits ${(data.limits || []).map((v) => v.toExponential(1)).join(" … ")}).`,
      "ok"
    );
  } catch (e) {
    setStatus(`Difference imaging request failed: ${e}`, "error");
  }
}

$("load").addEventListener("click", load);
$("view-orig").addEventListener("click", () => setViewMode("orig"));
$("view-diff").addEventListener("click", () => setViewMode("diff"));
$("preset").addEventListener("click", () => {
  $("ra").value = "133.795129";
  $("dec").value = "-7.245147";
});
$("play").addEventListener("click", () => {
  const list = activeFrames();
  if (!list.length) return;
  if (playing) { stop(); return; }
  if (current >= list.length - 1) showFrame(0);
  playing = true;
  $("play").textContent = "⏸ Pause";
  schedule();
});
$("slider").addEventListener("input", (e) => { stop(); showFrame(parseInt(e.target.value, 10)); });
$("speed").addEventListener("change", () => { if (playing) { stop(); playing = true; $("play").textContent = "⏸ Pause"; schedule(); } });
["ra", "dec", "fov"].forEach((id) => $(id).addEventListener("keydown", (e) => { if (e.key === "Enter") load(); }));
