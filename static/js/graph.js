var width = window.innerWidth * 0.75,
    height = 900;

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
        infoParagraph.innerHTML = data.id.replace("pm", "") + "<br>" + data.date_created + "<br>" + data.abstract;
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

    simulation.nodes(graph.nodes);
    simulation.force("link").links(graph.links);

    link = link
        .data(graph.links)
        .enter().append("line")
            .attr("class", function (d) {
                return "link " + d.type;
            });

    node = node
        .data(graph.nodes)
        .enter().append("circle")
            .attr("class", function(d) {
                return "node " + d.type;
            })
            .attr("r", 10)
            .call(d3.drag()
                .on("start", dragstarted)
                .on("drag", dragged)
                .on("end", dragended))
             .on("mouseover", function(d) {
                displayInfo(d);
             });
};
