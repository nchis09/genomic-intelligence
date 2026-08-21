// D3.js parallel coordinates for the nf-genomic-intelligence dashboard.
// Draws each inter-axis segment as a separate SVG line so it can be colored
// by the 0-100 score of the destination metric.

globalThis.parcoords = function(r2d3) {
  const data = r2d3.data;
  const samples = data.samples;
  const titles = data.titles;
  const darkMode = data.dark_mode === true;
  const nMetrics = titles.length;

  const margin = { top: 40, right: 100, bottom: 40, left: 70 };
  const width = r2d3.width - margin.left - margin.right;
  const height = r2d3.height - margin.top - margin.bottom;

  const textColor = darkMode ? "#e0e0e0" : "#333333";
  const gridColor = darkMode ? "#555555" : "#cccccc";

  r2d3.svg.selectAll("*").remove();

  const svg = r2d3.svg
    .attr("width", r2d3.width)
    .attr("height", r2d3.height)
    .append("g")
    .attr("transform", `translate(${margin.left},${margin.top})`);

  // Build per-metric y-scales from raw values
  const yScales = [];
  const metricDomains = [];
  for (let i = 0; i < nMetrics; i++) {
    const vals = samples.map(d => d.values[i]).filter(v => v !== null && !isNaN(v));
    const extent = d3.extent(vals);
    metricDomains.push(extent);
    yScales[i] = d3.scaleLinear()
      .domain(extent)
      .range([height, 0])
      .nice();
  }

  const xScale = d3.scalePoint()
    .domain(titles)
    .range([0, width])
    .padding(0.15);

  // Global 0-100 color scale (red -> yellow -> green -> cyan -> blue)
  const colorScale = d3.scaleSequential()
    .domain([0, 100])
    .interpolator(d3.interpolateRgbBasis([
      "#ff0000", "#ffff00", "#00ff00", "#00ffff", "#0000ff"
    ]));

  // Grid lines
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
    .attr("stroke-width", 1);

  // Axes
  const axisGroup = svg.append("g").attr("class", "axes");
  for (let i = 0; i < nMetrics; i++) {
    const x = xScale(titles[i]);
    const axis = d3.axisLeft(yScales[i]).ticks(5).tickSize(4);
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

  // Build segment list: one per sample per inter-axis pair
  const segments = [];
  samples.forEach(s => {
    for (let i = 0; i < nMetrics - 1; i++) {
      const x1 = xScale(titles[i]);
      const x2 = xScale(titles[i + 1]);
      const y1 = yScales[i](s.values[i]);
      const y2 = yScales[i + 1](s.values[i + 1]);
      const score = s.scores[i + 1]; // color by destination metric score
      segments.push({
        x1: x1,
        x2: x2,
        y1: y1,
        y2: y2,
        color: colorScale(score),
        sample: s.name
      });
    }
  });

  // Draw segments
  svg.selectAll(".segment")
    .data(segments)
    .enter()
    .append("line")
    .attr("class", "segment")
    .attr("x1", d => d.x1)
    .attr("x2", d => d.x2)
    .attr("y1", d => d.y1)
    .attr("y2", d => d.y2)
    .attr("stroke", d => d.color)
    .attr("stroke-width", 2)
    .attr("stroke-opacity", 0.85)
    .append("title")
    .text(d => d.sample);

  // Colorbar on the right
  const cbWidth = 16;
  const cbHeight = 150;
  const cbX = width + margin.right - 70;
  const cbY = (height - cbHeight) / 2;
  const nStops = 20;

  const defs = svg.append("defs");
  const gradient = defs.append("linearGradient")
    .attr("id", "parcoords-gradient")
    .attr("x1", "0%")
    .attr("y1", "100%")
    .attr("x2", "0%")
    .attr("y2", "0%");

  for (let i = 0; i <= nStops; i++) {
    const t = i / nStops;
    const score = t * 100;
    gradient.append("stop")
      .attr("offset", `${t * 100}%`)
      .attr("stop-color", colorScale(score));
  }

  const colorbar = svg.append("g").attr("class", "colorbar");

  colorbar.append("rect")
    .attr("x", cbX)
    .attr("y", cbY)
    .attr("width", cbWidth)
    .attr("height", cbHeight)
    .attr("fill", "url(#parcoords-gradient)")
    .attr("stroke", textColor)
    .attr("stroke-width", 0.5);

  // Colorbar title
  colorbar.append("text")
    .attr("x", cbX + cbWidth / 2)
    .attr("y", cbY - 10)
    .attr("text-anchor", "middle")
    .attr("fill", textColor)
    .style("font-size", "11px")
    .style("font-weight", "bold")
    .text("0-100 score");

  // Colorbar ticks
  const colorbarScale = d3.scaleLinear().domain([0, 100]).range([cbY + cbHeight, cbY]);
  const cbAxis = d3.axisRight(colorbarScale).ticks(5).tickSize(4);
  colorbar.append("g")
    .attr("transform", `translate(${cbX + cbWidth},0)`)
    .call(cbAxis)
    .selectAll("text")
    .attr("fill", textColor)
    .style("font-size", "10px");

  colorbar.selectAll("path, line")
    .attr("stroke", textColor);
};
