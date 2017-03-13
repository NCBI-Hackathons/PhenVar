var width = window.innerWidth * 0.75,
    height = window.innerHeight;

var svg = d3.select("svg")
    .attr("width", width)
    .attr("height", height)
    .call(d3.zoom().on("zoom", function () {
        svg.attr("transform", d3.event.transform)
    }))
    .append("g")
    .attr("width", width)
    .attr("height", height);

var simulation = d3.forceSimulation()
    .force("charge", d3.forceManyBody().strength(-150))
    .force("link", d3.forceLink().id(function(d) {return d.id;}).distance(100))
    .force("x", d3.forceX(width/2))
    .force("y", d3.forceY(height/2))
    .on("tick", ticked);

var link = svg.selectAll(".link"),
    node = svg.selectAll(".node");

var node_drag = d3.behavior.drag()
    .on("dragstart", dragstart)
    .on("drag", dragmove)
    .on("dragend", dragend);

function dragstarted(d) {
    if (!d3.event.active) simulation.alphaTarget(0.3).restart();
    d.fx = d.x;
    d.fy = d.y;
}

function displayInfo(data) {
    var infoParagraph = document.getElementById("info_paragraph");
    if (data.type == "article") {
        infoParagraph.innerHTML = data.id.replace("pm", "");
    }
    else {
        infoParagraph.innerHTML = data.id;
    }
}

function dragged(d) {
    d.fx = d3.event.x;
    d.fy = d3.event.y;
}

function dragended(d) {
    if (!d3.event.active) simulation.alphaTarget(0);
    d.fx = null;
    d.fy = null;
}

function ticked() {
    link.attr("x1", function(d) { return d.source.x; })
        .attr("y1", function(d) { return d.source.y; })
        .attr("x2", function(d) { return d.target.x; })
        .attr("y2", function(d) { return d.target.y; });
    node.attr("cx", function(d) { return d.x })
        .attr("cy", function(d) { return d.y });
}

function renderGraph(error, graph) {
    if (error) throw error;

    var maxWeight = Math.max.apply(Math, graph.nodes.map(function(node) {
        if (node.type == "article") {
            return node.weight;
        }
        else { return 0; }
    }));
    var maxSize = 24;
    var minSize = 6;
    var stepSize = (maxSize - minSize)/maxWeight

    console.log(maxWeight);
    simulation.nodes(graph.nodes);
    simulation.force("link").links(graph.links);

    link = link
        .data(graph.links)
        .enter().append("line")
            .attr("class", "link");

    node = node
        .data(graph.nodes)
        .enter().append("circle")
            .attr("class", "node")
            .attr("r", function(d) {  //10)
                if (d.type == "article") {
                    return d.weight * stepSize + minSize;
                }
                else {
                    return 6;
                }
            })
            .style("fill", function(d) {
                return {
                    "article": "blue",
                    "rsid": "red",
                }[d.type]
            })
            .call(d3.drag()
                .on("start", dragstarted)
                .on("drag", dragged)
                .on("end", dragended))
             .on("mouseover", function(d) {
                displayInfo(d);
             });
};
