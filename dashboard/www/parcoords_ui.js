// D3.js parallel coordinates for the nf-genomic-intelligence dashboard.
// Draws one smooth curved path per sample, colored per-sample.
// Plots raw metric values directly on a shared, static 0-100 y-axis.
// Includes a per-metric IQR shaded band and a checkbox filter sidebar
// that highlights samples even when curves overlap.

globalThis.drawParcoords = function(containerId, data) {
  const container = d3.select("#" + containerId);
  if (container.empty()) {
    console.error("parcoords container not found:", containerId);
    return;
  }

  container.selectAll("*").remove();
  container
    .style("display", "flex")
    .style("align-items", "flex-start")
    .style("gap", "16px");

  const samples = data.samples;
  const titles = data.titles;
  const darkMode = data.dark_mode === true;
  const nMetrics = titles.length;

  const rawWidth = Math.max(960, nMetrics * 220);
  const rawHeight = 520;
  const margin = { top: 40, right: 40, bottom: 40, left: 70 };
  const width = rawWidth - margin.left - margin.right;
  const height = rawHeight - margin.top - margin.bottom;

  const textColor = darkMode ? "#e0e0e0" : "#333333";
  const gridColor = darkMode ? "#555555" : "#cccccc";
  const bandFill = darkMode ? "rgba(180,180,200,0.14)" : "rgba(90,90,110,0.10)";

  const chartWrapper = container.append("div")
    .style("flex", "1 1 auto")
    .style("overflow-x", "auto");

  const svg = chartWrapper
    .append("svg")
    .attr("width", rawWidth)
    .attr("height", rawHeight)
    .attr("viewBox", `0 0 ${rawWidth} ${rawHeight}`)
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);

  // Shared static 0-100 y-scale for all metrics (raw values plotted directly)
  const sharedYScale = d3.scaleLinear()
    .domain([0, 100])
    .range([height, 0]);
  const yScales = titles.map(() => sharedYScale);

  const sampleByName = new Map(samples.map(s => [s.name, s]));
  const defaultColors = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"];

  function sampleColor(s, idx) {
    return s.color || defaultColors[idx % defaultColors.length];
  }

  function formatRaw(v) {
    return (v === null || typeof v === "undefined" || isNaN(v)) ? "NA" : d3.format(".4g")(v);
  }

  const xScale = d3.scalePoint()
    .domain(titles)
    .range([0, width])
    .padding(0.15);

  // Grid lines (faint, recede behind the data)
  svg.selectAll(".grid-line")
    .data(d3.range(nMetrics))
    .enter()
    .append("line")
    .attr("class", "grid-line")
    .attr("x1", d => xScale(titles[d]))
    .attr("x2", d => xScale(titles[d]))
    .attr("y1", 0)
    .attr("y2", height)
    .attr("stroke", gridColor)
    .attr("stroke-width", 1)
    .attr("stroke-opacity", 0.35);

  // Shaded bands per metric (help spot outliers), drawn behind curves:
  // a wider mean +/- 2 SD "normal range" halo, with the tighter IQR band on top.
  function quantileOf(vals, p) {
    const sorted = vals.filter(v => v !== null && !isNaN(v)).slice().sort((a, b) => a - b);
    if (sorted.length === 0) return 0;
    return d3.quantileSorted(sorted, p);
  }

  function meanAndStdDev(vals) {
    const clean = vals.filter(v => v !== null && !isNaN(v));
    if (clean.length === 0) return { mean: 0, sd: 0 };
    const mean = d3.mean(clean);
    const sd = clean.length > 1 ? d3.deviation(clean) : 0;
    return { mean: mean, sd: sd || 0 };
  }

  const bandData = titles.map((t, i) => {
    const vals = samples.map(s => s.values[i]);
    const { mean, sd } = meanAndStdDev(vals);
    return {
      q1: quantileOf(vals, 0.25),
      q3: quantileOf(vals, 0.75),
      sdLow: Math.max(0, mean - 2 * sd),
      sdHigh: Math.min(100, mean + 2 * sd)
    };
  });

  const iqrAreaGen = d3.area()
    .x((d, i) => xScale(titles[i]))
    .y0(d => sharedYScale(d.q1))
    .y1(d => sharedYScale(d.q3))
    .curve(d3.curveCatmullRom.alpha(0.5));

  const sdAreaGen = d3.area()
    .x((d, i) => xScale(titles[i]))
    .y0(d => sharedYScale(d.sdLow))
    .y1(d => sharedYScale(d.sdHigh))
    .curve(d3.curveCatmullRom.alpha(0.5));

  const sdBandFill = darkMode ? "rgba(150,150,180,0.07)" : "rgba(90,90,110,0.05)";

  svg.append("path")
    .datum(bandData)
    .attr("class", "sd-band")
    .attr("d", sdAreaGen)
    .attr("fill", sdBandFill)
    .attr("stroke", "none");

  svg.append("path")
    .datum(bandData)
    .attr("class", "iqr-band")
    .attr("d", iqrAreaGen)
    .attr("fill", bandFill)
    .attr("stroke", "none");

  // Axes -- tick number labels are only shown on the leftmost axis; the
  // remaining axes show tick marks/lines without repeated numbers.
  const axisGroup = svg.append("g").attr("class", "axes");
  for (let i = 0; i < nMetrics; i++) {
    const x = xScale(titles[i]);
    const axis = d3.axisLeft(yScales[i]).ticks(5).tickSize(4);
    if (i > 0) axis.tickFormat("");
    const g = axisGroup.append("g")
      .attr("transform", `translate(${x},0)`)
      .call(axis);

    g.selectAll("text")
      .attr("fill", textColor)
      .style("font-size", "11px");
    g.selectAll("path, line")
      .attr("stroke", textColor);

    // Axis title
    axisGroup.append("text")
      .attr("x", x)
      .attr("y", -20)
      .attr("text-anchor", "middle")
      .attr("fill", textColor)
      .style("font-size", "12px")
      .style("font-weight", "bold")
      .text(titles[i]);
  }

  // Build one smooth curved path per sample (raw values)
  const lineGen = d3.line()
    .x((d, i) => xScale(titles[i]))
    .y((d, i) => yScales[i](d))
    .curve(d3.curveCatmullRom.alpha(0.5));

  const paths = samples.map((s, idx) => ({
    sample: s.name,
    color: sampleColor(s, idx),
    d: lineGen(s.values)
  }));

  // Popup elements (created once per chart, hidden by default)
  const popup = container.append("div")
    .attr("class", "parcoords-popup")
    .style("position", "fixed")
    .style("display", "none")
    .style("z-index", "9999")
    .style("padding", "12px 14px")
    .style("border-radius", "8px")
    .style("border", `1px solid ${gridColor}`)
    .style("background", darkMode ? "rgba(30,30,30,0.97)" : "rgba(255,255,255,0.97)")
    .style("color", textColor)
    .style("font-family", "sans-serif")
    .style("font-size", "12px")
    .style("box-shadow", "0 4px 14px rgba(0,0,0,0.3)")
    .style("max-width", "420px");
  popup.on("click", event => event.stopPropagation());

  let pinnedSample = null;
  const highlightedSamples = new Set();

  function setCheckboxState() {
    sidebar.selectAll("input[type=checkbox]").each(function() {
      this.checked = highlightedSamples.has(this.getAttribute("data-sample"));
    });
  }

  function applyHighlight() {
    if (highlightedSamples.size === 0) {
      svg.selectAll(".sample-path")
        .attr("stroke-opacity", 0.85)
        .attr("stroke-width", 2);
    } else {
      svg.selectAll(".sample-path")
        .attr("stroke-opacity", d => highlightedSamples.has(d.sample) ? 1 : 0.12)
        .attr("stroke-width", d => highlightedSamples.has(d.sample) ? 3 : 1.5);
    }
  }

  function metricBox(indices) {
    const rows = indices.map(i => `<div style="margin-bottom:6px">` +
      `<div style="font-weight:600">${titles[i]}</div>` +
      `<div>value: ${formatRaw(sampleByName.get(pinnedSample).values[i])}</div>` +
      `</div>`).join("");
    return `<div style="flex:1; min-width:150px; padding:8px; border:1px solid ${gridColor}; border-radius:6px;">${rows}</div>`;
  }

  function showPopup(sample, event) {
    pinnedSample = sample;
    const s = sampleByName.get(sample);
    if (!s) return;
    const half = Math.ceil(nMetrics / 2);
    const idxAll = d3.range(nMetrics);
    const box1 = metricBox(idxAll.slice(0, half));
    const box2 = metricBox(idxAll.slice(half));
    const refLine = s.closest_reference
      ? `<div style="font-size:12.5px; font-weight:500; color:${textColor}; margin-top:2px;">Closest reference: ${s.closest_reference}</div>`
      : "";
    popup.html(
      `<div style="display:flex; justify-content:space-between; align-items:flex-start; margin-bottom:8px;">` +
        `<div><div style="font-weight:bold;">${s.name}</div>${refLine}</div>` +
        `<div class="parcoords-popup-close" style="cursor:pointer; font-weight:bold; padding:0 4px;">&times;</div>` +
      `</div>` +
      `<div style="display:flex; gap:10px;">${box1}${box2}</div>`
    );
    popup.select(".parcoords-popup-close").on("click", function(evt) {
      evt.stopPropagation();
      closePopup();
    });
    positionPopup(event);
    popup.style("display", "block");
    highlightedSamples.clear();
    highlightedSamples.add(sample);
    setCheckboxState();
    applyHighlight();
  }

  function positionPopup(event) {
    const vw = window.innerWidth;
    const vh = window.innerHeight;
    let left = event.clientX + 14;
    let top = event.clientY + 14;
    const popupWidth = 380;
    const popupHeight = 220;
    if (left + popupWidth > vw) left = event.clientX - popupWidth - 14;
    if (top + popupHeight > vh) top = event.clientY - popupHeight - 14;
    popup.style("left", left + "px").style("top", top + "px");
  }

  function closePopup() {
    pinnedSample = null;
    popup.style("display", "none");
    highlightedSamples.clear();
    setCheckboxState();
    applyHighlight();
  }

  // Draw one curved path per sample
  svg.selectAll(".sample-path")
    .data(paths)
    .enter()
    .append("path")
    .attr("class", "sample-path")
    .attr("d", d => d.d)
    .attr("fill", "none")
    .attr("stroke", d => d.color)
    .attr("stroke-width", 2)
    .attr("stroke-opacity", 0.85)
    .attr("stroke-linecap", "round")
    .attr("stroke-linejoin", "round")
    .style("cursor", "pointer")
    .on("click", function(event, d) {
      event.stopPropagation();
      if (pinnedSample === d.sample) {
        closePopup();
      } else {
        pinnedSample = d.sample;
        showPopup(d.sample, event);
      }
    })
    .append("title")
    .text(d => d.sample);

  // Click on chart background (but not the popup or sidebar) closes the popup
  container.on("click", function() {
    closePopup();
  });

  // Filter sidebar: checkboxes highlight one or more samples even when
  // their curves overlap on the chart.
  const sidebar = container.append("div")
    .attr("class", "parcoords-filter-sidebar")
    .style("flex", "0 0 170px")
    .style("padding", "10px 12px")
    .style("border-radius", "8px")
    .style("border", `1px solid ${gridColor}`)
    .style("background", darkMode ? "rgba(255,255,255,0.03)" : "rgba(0,0,0,0.02)")
    .style("font-family", "sans-serif")
    .style("font-size", "12px")
    .style("color", textColor);
  sidebar.on("click", event => event.stopPropagation());

  sidebar.append("div")
    .style("font-weight", "bold")
    .style("margin-bottom", "8px")
    .text("Highlight samples");

  sidebar.append("input")
    .attr("type", "search")
    .attr("placeholder", "Search samples...")
    .style("width", "100%")
    .style("box-sizing", "border-box")
    .style("padding", "5px 8px")
    .style("margin-bottom", "8px")
    .style("border-radius", "4px")
    .style("border", `1px solid ${gridColor}`)
    .style("background", darkMode ? "rgba(0,0,0,0.2)" : "#ffffff")
    .style("color", textColor)
    .style("font-size", "12px")
    .on("input", function() {
      const q = this.value.trim().toLowerCase();
      rowsGroup.selectAll(".parcoords-filter-row")
        .style("display", d => d.name.toLowerCase().includes(q) ? "flex" : "none");
    });

  const rowsGroup = sidebar.append("div")
    .style("max-height", "320px")
    .style("overflow-y", "auto");

  samples.forEach((s, idx) => {
    const color = sampleColor(s, idx);
    const row = rowsGroup.append("label")
      .datum(s)
      .attr("class", "parcoords-filter-row")
      .style("display", "flex")
      .style("align-items", "center")
      .style("gap", "6px")
      .style("margin-bottom", "6px")
      .style("cursor", "pointer");

    row.append("input")
      .attr("type", "checkbox")
      .attr("data-sample", s.name)
      .style("margin", "0")
      .on("change", function() {
        if (this.checked) highlightedSamples.add(s.name);
        else highlightedSamples.delete(s.name);
        pinnedSample = null;
        popup.style("display", "none");
        applyHighlight();
      });

    row.append("span")
      .style("width", "10px")
      .style("height", "10px")
      .style("border-radius", "2px")
      .style("flex", "0 0 auto")
      .style("background", color);

    row.append("span").text(s.name);
  });
};
