import random as rand
import tkinter as tk
from math import exp, log
from tkinter import ttk

import matplotlib.pyplot
from PIL import Image, ImageTk
from graphviz import Graph

seed = 4587955
alpha = 0.5


class MainWindow:
    root = tk.Tk()
    img_label = None
    img_count = 0
    img_max = -1
    node_txt, edge_txt, p0_txt, err_lbl = None, None, None, None
    nodes = None
    edges = None
    adj_lists = None
    p0 = None

    def __init__(self, nodes: int, edges: int, length: int, adj_lists: [], p0: int = 0):
        self.nodes = nodes
        self.edges = edges
        self.img_max = length
        self.adj_lists = adj_lists
        self.p0 = p0
        title = "Tutorial 3: Techniques for Modelling Contact Networks in a Pandemic"
        self.root.title(title)
        frame = tk.Frame(self.root)

        image = Image.open("graph" + str(self.img_count) + ".png")
        wi, hi = image.size
        wi, hi = size_check(wi, hi)
        image = image.resize((wi, hi), Image.ANTIALIAS)
        photo = ImageTk.PhotoImage(image)

        self.img_label = tk.Label(frame, image=photo, bg='green', anchor="center")
        self.img_label.image = photo
        self.img_label.grid(row=0, columnspan=2, padx=2, pady=2)

        prev_btn = tk.Button(frame, text="Previous Day", anchor="center", command=self.prev_img)
        prev_btn.grid(row=1, column=0, padx=2, pady=2, sticky="E")

        prev_btn = tk.Button(frame, text="Next Day", anchor="center", command=self.next_img)
        prev_btn.grid(row=1, column=1, padx=2, pady=2, sticky="W")

        node_lbl = tk.Label(frame, text="Number of Nodes: ", anchor="center")
        node_lbl.grid(row=2, column=0, padx=2, pady=2, sticky="E")

        self.node_txt = tk.Text(frame, height=1, width=5)
        self.node_txt.grid(row=2, column=1, padx=2, pady=2, sticky="W")
        self.node_txt.insert('end', "64")

        edge_lbl = tk.Label(frame, text="Number of Edges: ", anchor="center")
        edge_lbl.grid(row=3, column=0, padx=2, pady=2, sticky="E")

        self.edge_txt = tk.Text(frame, height=1, width=5, wrap='none')
        self.edge_txt.grid(row=3, column=1, padx=2, pady=2, sticky="W")
        self.edge_txt.insert('end', "128")

        p0_lbl = tk.Label(frame, text="Patient Zero: ", anchor="center")
        p0_lbl.grid(row=4, column=0, padx=2, pady=2, sticky="E")

        self.p0_txt = tk.Text(frame, height=1, width=5, wrap='none')
        self.p0_txt.grid(row=4, column=1, padx=2, pady=2, sticky="W")
        self.p0_txt.insert('end', "0")

        new_p0_btn = tk.Button(frame, text="Change Patient Zero", anchor="center", command=self.update_patient0)
        new_p0_btn.grid(row=5, column=0, padx=2, pady=2, sticky="E")

        gen_btn = tk.Button(frame, text="Generate New Random Network", anchor="center", command=self.gen_network)
        gen_btn.grid(row=5, column=1, padx=2, pady=2, sticky="W")

        PLC_btn = tk.Button(frame, text="Generate New PLC Network", anchor="center", command=self.gen_PLC)
        PLC_btn.grid(row=6, column=0, padx=2, pady=2, sticky="E")

        ER_btn = tk.Button(frame, text="Generate New ER Network", anchor="center", command=self.gen_ER)
        ER_btn.grid(row=6, column=1, padx=2, pady=2, sticky="W")

        WS_btn = tk.Button(frame, text="Generate New WS Network", anchor="center", command=self.gen_WS)
        WS_btn.grid(row=7, column=0, padx=2, pady=2, sticky="E")

        sim_btn = tk.Button(frame, text="Simulate 50 Epidemics", anchor="center", command=self.simulate_epis)
        sim_btn.grid(row=7, column=1, padx=2, pady=2, sticky="W")

        ttk.Separator(frame, orient='horizontal').grid(row=8, columnspan=2, sticky="EW")

        self.err_lbl = tk.Label(frame, text="", anchor="center")
        self.err_lbl.grid(row=9, columnspan=2, padx=2, pady=2, sticky="EW")

        frame.grid(row=0, column=0, sticky="NSEW")
        self.root.mainloop()
        pass

    def prev_img(self):
        if self.img_count > 0:
            self.img_count -= 1
            pass
        self.change_img()
        self.err_lbl.config(text="")
        pass

    def next_img(self):
        if self.img_count < self.img_max - 1:
            self.img_count += 1
            pass
        self.change_img()
        self.err_lbl.config(text="")
        pass

    def change_img(self):
        image = Image.open("graph" + str(self.img_count) + ".png")
        wi, hi = image.size
        wi, hi = size_check(wi, hi)
        image = image.resize((wi, hi), Image.ANTIALIAS)
        photo = ImageTk.PhotoImage(image)
        self.img_label.configure(image=photo)
        self.img_label.image = photo
        self.root.update_idletasks()

    def gen_network(self):
        if self.txt_check():
            self.nodes = int(self.node_txt.get("1.0", 'end-1c'))
            self.edges = int(self.edge_txt.get("1.0", 'end-1c'))
            self.p0 = int(self.p0_txt.get("1.0", 'end-1c'))
            self.adj_lists, nodes = make_network("edgelists.txt", self.nodes, self.edges)
            edg_list = get_edge_list("edgelists.txt")
            epi_log = run_epi(self.adj_lists, self.nodes, self.p0)
            make_graphs(edg_list, epi_log)
            self.img_count = 0
            self.img_max = int(len(epi_log))
            self.change_img()
            self.err_lbl.config(text="")
            pass
        pass

    def gen_PLC(self):
        if self.txt_check():
            self.nodes = int(self.node_txt.get("1.0", 'end-1c'))
            self.edges = int(self.edge_txt.get("1.0", 'end-1c'))
            self.p0 = int(self.p0_txt.get("1.0", 'end-1c'))
            self.adj_lists, nodes = powerlaw_cluster("edgelists.txt", self.nodes, self.edges, 0.2)
            edg_list = get_edge_list("edgelists.txt")
            epi_log = run_epi(self.adj_lists, self.nodes, self.p0)
            make_graphs(edg_list, epi_log)
            self.img_count = 0
            self.img_max = int(len(epi_log))
            self.change_img()
            self.err_lbl.config(text="")
            pass
        pass

    def gen_ER(self):
        if self.txt_check():
            self.nodes = int(self.node_txt.get("1.0", 'end-1c'))
            self.edges = int(self.edge_txt.get("1.0", 'end-1c'))
            self.p0 = int(self.p0_txt.get("1.0", 'end-1c'))
            self.adj_lists, nodes = erdos_renyi("edgelists.txt", self.nodes, 0.05)
            edg_list = get_edge_list("edgelists.txt")
            epi_log = run_epi(self.adj_lists, self.nodes, self.p0)
            make_graphs(edg_list, epi_log)
            self.img_count = 0
            self.img_max = int(len(epi_log))
            self.change_img()
            self.err_lbl.config(text="")
            pass
        pass

    def gen_WS(self):
        if self.txt_check():
            self.nodes = int(self.node_txt.get("1.0", 'end-1c'))
            self.edges = int(self.edge_txt.get("1.0", 'end-1c'))
            self.p0 = int(self.p0_txt.get("1.0", 'end-1c'))
            self.adj_lists, nodes = watts_stogatz("edgelists.txt", self.nodes, 4, 0.2)
            edg_list = get_edge_list("edgelists.txt")
            epi_log = run_epi(self.adj_lists, self.nodes, self.p0)
            make_graphs(edg_list, epi_log)
            self.img_count = 0
            self.img_max = int(len(epi_log))
            self.change_img()
            self.err_lbl.config(text="")
            pass
        pass

    def update_patient0(self):
        if self.txt_check():
            self.p0 = int(self.p0_txt.get("1.0", 'end-1c'))
            edg_list = get_edge_list("edgelists.txt")
            epi_log = run_epi(self.adj_lists, self.nodes, self.p0)
            make_graphs(edg_list, epi_log)
            self.img_count = 0
            self.img_max = int(len(epi_log))
            self.change_img()
            self.err_lbl.config(text="")
            pass
        pass

    def txt_check(self):
        nodes = int(self.node_txt.get("1.0", 'end-1c'))
        edges = int(self.edge_txt.get("1.0", 'end-1c'))
        p0 = int(self.p0_txt.get("1.0", 'end-1c'))
        if nodes > edges:
            self.err_lbl.config(text="Nodes must be less than or equal to the number of edges.")
            self.root.update_idletasks()
            return False
        elif p0 >= nodes:
            self.err_lbl.config(text="Patient zero must be less than the number of nodes.")
            self.root.update_idletasks()
            return False
        return True

    def simulate_epis(self):
        run_epis(self.adj_lists, self.nodes, self.p0)
        image = Image.open("epi_profile.png")
        wi, hi = image.size
        wi, hi = size_check(wi, hi)
        image = image.resize((wi, hi), Image.ANTIALIAS)
        photo = ImageTk.PhotoImage(image)
        self.img_label.configure(image=photo)
        self.img_label.image = photo
        self.img_count = 0
        self.root.update_idletasks()
        pass


def size_check(wi, hi):
    lrg = max(wi, hi)
    nlrg = 700
    ratio = nlrg / lrg
    return int(wi * ratio), int(hi * ratio)


def make_network(filename: str, nodes: int, edges: int):
    adj_lists = [[] for _ in range(nodes)]
    edg_count = 0
    for n in range(nodes):
        while True:
            to = rand.randint(0, nodes - 1)
            if to != n:
                break
                pass
            pass
        # print(str(n) + '\t' + str(to))
        adj_lists[n].append(to)
        adj_lists[to].append(n)
        edg_count += 1
        pass

    edg_to_add = edges - edg_count
    for _ in range(edg_to_add):
        while True:
            fr = rand.randint(0, nodes - 1)
            to = rand.randint(0, nodes - 1)
            if fr != to:
                if fr not in adj_lists[to]:
                    break
                    pass
                pass
            pass
        adj_lists[fr].append(to)
        adj_lists[to].append(fr)
        edg_count += 1
        pass

    with open(filename, "w") as f:
        f.write("Nodes: " + str(nodes) + '\t')
        f.write("Edges: " + str(edg_count) + '\n')
        for li in adj_lists:
            for n in li:
                f.write(str(n) + '\t')
                pass
            f.write('\n')
            pass
        pass
    return adj_lists, nodes


def erdos_renyi(filename: str, nodes: int, edge_probability: float):
    adj_lists = [[] for _ in range(nodes)]
    edg_count = 0
    for n in range(nodes):
        for m in range(nodes):
            if n == m:
                continue
            if rand.random() <= (1 - edge_probability):
                continue
            adj_lists[n].append(m)
            adj_lists[m].append(n)
            edg_count += 1
            pass
        pass

    with open(filename, "w") as f:
        f.write("Nodes: " + str(nodes) + '\t')
        f.write("Edges: " + str(edg_count) + '\n')
        for li in adj_lists:
            for n in li:
                f.write(str(n) + '\t')
                pass
            f.write('\n')
            pass
        pass
    return adj_lists, nodes


def watts_stogatz(filename: str, nodes: int, k: int, beta: float):
    adj_lists = [[] for _ in range(nodes)]
    edg_count = 0
    for n in range(nodes):
        for m in range(int(k / 2)):
            if m == int(k / 2) - 1:
                continue
            to = (m + 1 + n) % nodes
            adj_lists[n].append(to)
            adj_lists[to].append(n)
            edg_count += 1
            pass
        pass

    for n in range(nodes):
        if rand.random() > beta:
            to = int(k / 2 + n) % nodes
            adj_lists[n].append(to)
            adj_lists[to].append(n)
            edg_count += 1
            continue
        while True:
            fr = rand.randint(0, nodes - 1)
            if fr != n:
                if fr not in adj_lists[n]:
                    break
                    pass
                pass
            pass
        adj_lists[fr].append(n)
        adj_lists[n].append(fr)
        edg_count += 1
        pass

    with open(filename, "w") as f:
        f.write("Nodes: " + str(nodes) + '\t')
        f.write("Edges: " + str(edg_count) + '\n')
        for li in adj_lists:
            for n in li:
                f.write(str(n) + '\t')
                pass
            f.write('\n')
            pass
        pass
    return adj_lists, nodes


def powerlaw_cluster(filename: str, nodes: int, edges: int, beta: float):
    adj_lists = [[] for _ in range(nodes)]
    edg_count = 0
    to = 0
    for _ in range(edges):
        while True:
            fr = rand.randint(0, nodes - 1)
            to = rand.randint(0, nodes - 1)
            if fr != to:
                if fr not in adj_lists[to]:
                    break
                    pass
                pass
            pass
        adj_lists[fr].append(to)
        adj_lists[to].append(fr)
        edg_count += 1
        if edg_count == 1:
            pass
        fr2 = 0
        if rand.random() <= beta:
            if len(adj_lists[fr]) > 0:
                fr2 = adj_lists[fr][rand.randint(0, len(adj_lists[fr]) - 1)]
                pass
            pass
        if fr2 not in adj_lists[to]:
            if fr2 != to:
                adj_lists[fr2].append(to)
                adj_lists[to].append(fr2)
                edg_count += 1
                pass
            pass
        pass

    with open(filename, "w") as f:
        f.write("Nodes: " + str(nodes) + '\t')
        f.write("Edges: " + str(edg_count) + '\n')
        for li in adj_lists:
            for n in li:
                f.write(str(n) + '\t')
                pass
            f.write('\n')
            pass
        pass
    return adj_lists, nodes


def make_graph(inp: [], out: str, inf: [], rem: []):
    g = Graph(engine='sfdp')
    g.attr(size="6,6")
    g.graph_attr.update(dpi='600', overlap='false')
    g.node_attr.update(shape='circle', style='filled', fontsize='12', fixedsize='true')
    g.edge_attr.update()

    for n in range(inp[0]):
        if n in inf:
            g.node(str(n), label=str(n), fillcolor='red')
            pass
        elif n in rem:
            g.node(str(n), label=str(n), fillcolor='green')
            pass
        else:
            g.node(str(n), label=str(n), fillcolor='white')
            pass
        pass

    for e in inp[1:]:
        if e[0] >= e[1]:
            g.edge(str(e[0]), str(e[1]))
            pass
        pass

    g.render(filename=out, cleanup=True, format='png')
    pass


def get_edge_list(inp: str) -> []:
    with open(inp, "r") as f:
        first_line = f.readline()
        nodes = int(first_line.rstrip().split('\t')[0].split(' ')[1])
        edg_list = [nodes]
        lines = f.readlines()
        for fr, line in enumerate(lines):
            line = line.rstrip()
            line = line.split('\t')
            if len(line) > 0:
                for to in line:
                    if to != '':
                        edg_list.append([fr, int(to)])
                        pass
                    pass
                pass
            pass
    return edg_list


def make_graphs(inp: [], epi_log: []):
    removed = []
    for num, li in enumerate(epi_log):
        make_graph(inp, "graph" + str(num), li, removed)
        for n in li:
            removed.append(n)
            pass
        pass
    pass


def infected(sick: int):
    beta = 1 - exp(sick * log(1 - alpha))
    if rand.random() < beta:
        return True
    return False


def run_epis(adj_lists, nodes, p0: int = 0):
    epidemics = 50
    epi_logs = []
    lengths = []
    sums = [0 for _ in range(nodes)]
    counts = [0 for _ in range(nodes)]
    for _ in range(epidemics):
        epi_log = run_epi(adj_lists, nodes, p0)
        epi_log = [len(n) for n in epi_log]
        epi_logs.append(epi_log)
        lengths.append(len(epi_log))
        pass

    for ln, el in enumerate(epi_logs):
        for day in range(lengths[ln]):
            sums[day] += el[day]
            counts[day] += 1
            pass
        pass

    avg = []
    avg_all = []
    for day, s in enumerate(sums):
        if counts[day] > 0:
            avg.append(s / counts[day])
            avg_all.append(s / epidemics)
        pass

    max_len = max(lengths)
    for el in epi_logs:
        s_len = len(el)
        for _ in range(max_len - s_len):
            el.append(0)
    x = [n for n in range(max_len)]
    x_lbls = [str(n) for n in range(max_len)]

    fig = matplotlib.pyplot.figure()
    fig.set_dpi(400)
    fig.set_figheight(4)
    plot = fig.add_subplot(111)

    for el in epi_logs:
        plot.plot(x, el, linewidth=1, alpha=0.3, color='gray')
        pass

    plot.plot(x, avg, label="Average of Running")
    plot.plot(x, avg_all, label="Average of All")
    fig.suptitle("Epidemic Profiles for 50 Epidemics")
    plot.set_ylabel("Newely Infected Individuals")
    plot.set_xlabel("Day")
    plot.set_xticks(x)
    plot.set_xticklabels(x_lbls)
    plot.legend()
    fig.tight_layout()
    fig.savefig("epi_profile.png")


def run_epi(adj_lists: [], nodes: int, p0: int = 0):
    n_state = [0 for _ in range(nodes)]  # susceptible
    n_state[p0] = 1
    epi_log = [[p0]]
    num_infected = 1
    ttl_infected = 0

    length = 0
    while num_infected > 0:
        inf_neighbours = [0 for _ in range(nodes)]

        for n in range(nodes):
            if n_state[n] == 1:
                for nei in adj_lists[n]:
                    inf_neighbours[nei] += 1
                    pass
                pass
            pass

        for n in range(nodes):
            if n_state[n] == 0 and inf_neighbours[n] > 0:
                if infected(inf_neighbours[n]):
                    n_state[n] = 3
                    pass
                pass
            pass

        ttl_infected += num_infected
        num_infected = 0
        new_inf = []
        for n in range(nodes):
            if n_state[n] == 1:  # infected -> removed
                n_state[n] = 2
                pass
            elif n_state[n] == 3:
                n_state[n] = 1
                num_infected += 1
                new_inf.append(n)
                pass
            pass
        epi_log.append(new_inf)
        length += 1
        pass
    return epi_log


def main():
    rand.seed(seed)
    nodes = 64
    edges = 128
    adj_lists, nodes = make_network("edgelists.txt", nodes, edges)
    # adj_lists, nodes = powerlaw_cluster("edgelists.txt", nodes, edges, 0.2)
    # adj_lists, nodes = watts_stogatz("edgelists.txt", nodes, 4, 0.2)
    # adj_lists, nodes = erdos_renyi("edgelists.txt", nodes, 0.05)
    edg_list = get_edge_list("edgelists.txt")
    epi_log = run_epi(adj_lists, nodes)
    make_graphs(edg_list, epi_log)
    MainWindow(nodes, edges, len(epi_log), adj_lists)
    pass


main()
