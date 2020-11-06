#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <numeric>
using namespace std;

#define LOG(...) handleLog(#__VA_ARGS__, __VA_ARGS__)
#define LL(...)  std::cout << (__VA_ARGS__) << "\n"
#define LN       std::cout << "\n";

// кол-во маркеров для каждого атома
int markerCount = 3;
    
template <typename T>
void handleLogImpl(const T& value) {
    std::cout << value;
}

template <typename Head, typename ... Tail>
void handleLogImpl(const Head& head, const Tail&... tail) {
    std::cout << head << ", ";
    handleLogImpl(tail...);
}

template <typename ... Args>
void handleLog(const char* macroArgsAsStr, const Args&... args) {
    std::cout << macroArgsAsStr << " = ";
    handleLogImpl(args...);
    std::cout << "\n";
}

template <typename T>
ostream& operator<<(ostream& out, const vector<T>& vector) {
    out << "[";
    for (unsigned i = 0; i < vector.size(); i++) {
        out << vector[i];
        if (i != vector.size() - 1) {
            out << ", ";
        }
    }
    out << "]";
    return out;
}

template <typename Key, typename T>
ostream& operator<<(ostream& out, const map<Key, T>& map) {
    out << "{";
    unsigned i = 0;
    for (auto& it : map) {
        out << it.first << ": " << it.second;
        if (i++ != map.size() - 1) {
            out << ", ";
        }
    }
    out << "}";
    return out;
}

struct Link {
    int fst, snd;
    int type;
};

ostream& operator<<(ostream& out, Link l) {
    out << "fst=" << l.fst << ", snd=" << l.snd << ", type=" << l.type;
    return out;
}

struct Atom {
    int index;
    string name;
    int hydrogenCount;  // кол-во водорода (в mol файлах всегда 0)
    float x, y, z = 0;
};

ostream& operator<<(ostream& out, Atom a) {
    out << "index=" << a.index << ", name=" << a.name << ", hydrogenCount=" << a.hydrogenCount << ", ";
    out << "x,y,z={" << a.x << "," << a.y << "," << a.z << "}";
    return out;
}

struct Graph {
    int vertexCount;
    vector<vector<int>> table;
    
    Graph(int _vertexCount, const vector<vector<int>> _table) {
        vertexCount = _vertexCount;
        table = _table;
    }
    
    // возвращает соседей vertex'a
    vector<int> neighbors(int vertex) {
        vector<int> ret;
        for (int i = 0; i < vertexCount; i++) {
            if (table[vertex][i] != 0) {
                ret.push_back(i);
            }
        }
        return ret;
    }
    
    // ищет все подграфы vertex'a длиной distance (только для distance <= 3)
    vector<vector<int>> subgraph(int vertex, int distance) {
        if (distance == 1) {
            return {{vertex}};
        }
        
        else if (distance == 2) {
            vector<vector<int>> ret;
            for (auto& neighbor : neighbors(vertex)) {
                ret.push_back({vertex, neighbor});
            }
            return ret;
        }
        
        else if (distance == 3) {
            vector<vector<int>> ret;
            for (auto& neighbor1 : neighbors(vertex)) {
                for (auto& neighbor2 : neighbors(neighbor1)) {
                    if (neighbor2 != vertex) {
                        ret.push_back({vertex, neighbor1, neighbor2});
                    }
                }
            }
            return ret;
        }
        
        else {
            cout << "Graph::subgraph: distance > 3 (" << distance << " > 3)";
            exit(-1);
        }
        
        return {};
    }
    
    // ищет ВСЕ подграфы длиной distance
    vector<vector<int>> allSubgraphs(int distance) {
        vector<vector<int>> ret;
        for (int i = 0; i < vertexCount; i++) {
            vector<vector<int>> v = subgraph(i, distance);
            ret.insert(ret.end(), v.begin(), v.end());
        }
        for (int i = 0; i < int(ret.size()) - 1; i++) {
            vector<int> reversed = ret[i];
            reverse(reversed.begin(), reversed.end());
            for (int k = i + 1; k < int(ret.size()); k++) {
                if (ret[k] == reversed) {
                    ret.erase(ret.begin() + k);
                }
            }
        }
        return ret;
    }
    
    vector<int> reachableVertices(int i) {
        vector<int> visited(vertexCount);
        fill(visited.begin(), visited.end(), 0);
        wave(i, visited);
        return visited;
    }
    
    void wave(int i, vector<int>& visited) {
        if (visited[i] == 0) {
            visited[i]++;
            for (int k = 0; k < vertexCount; k++) {
                if (table[i][k] != 0) {
                    wave(k, visited);
                }
            }
        }
    }
};

ostream& operator<<(ostream& out, Graph g) {
    out << "\n";
    for (auto it : g.table) {
        out << it << ", " << endl;
    }
    return out;
}

struct Molecule {
    string name;
    int atomCount;
    vector<Atom> atoms;
    vector<Link> links;
    
    // создаёт граф, заполненный типами связей
    Graph createGraph() const {
        vector<vector<int>> table(atomCount);
        for (auto& v : table) {
            v.resize(atomCount);
            fill(v.begin(), v.end(), 0);
        }
        for (auto& l : links) {
            table[l.fst][l.snd] = l.type;
            table[l.snd][l.fst] = l.type;
        }
        return Graph(atomCount, table);
    }
    
    // находит связь между атомами a1, a2
    Link findLinkBetween(const Atom& a1, const Atom& a2) {
        int a = a1.index;
        int b = a2.index;
        for (const auto& link : links) {
            if ((link.fst == a && link.snd == b) || (link.fst == b && link.snd == a)) {
                return link;
            }
        }
        cout << "Molecule::findLinkBetween: can\'t find a link between " << a1 << " & " << a2 << endl;
        exit(-1);
    }
    
    // возвращает все связи атома под индексом atomIndex
    vector<Link> atomLinks(int atomIndex) {
        vector<Link> ret;
        for (auto it : links) {
            if (it.fst == atomIndex || it.snd == atomIndex) {
                ret.push_back(it);
            }
        }
        return ret;
    }
    
    // кол-во связей у атома с индексом atomIndex
    int edgeCount(int atomIndex) {
        return atomLinks(atomIndex).size();
    }
    
    // маркируем связь атома (одинарная, двойная и т.п.) (определяем буковку 's', 'd', ...)
    char markAtomLink(int atomIndex) {
        vector<Link> links = atomLinks(atomIndex);
        bool types[5] = {true, false, false, false, false};
        char syms[5] = {'s', 'd', 't', 'w', 'a'};
        for (auto it : links) {
            if (it.type != 1) types[0] = false;
            if (it.type == 2) types[1] = true;
            if (it.type == 3) types[2] = true;
            if (it.type == 4) types[3] = true;
            if (it.type == 5) types[4] = true;
        }
        for (int i = 0; i < 5; i++) {
            if (types[i]) {
                return syms[i];
            }
        }
        return '?';
    }
    
    string markAtom(int atomIndex) {
        string ret;
        string atomName = atoms[atomIndex].name;
        
        if (atomName.size() == 2) ret += atomName;
        else ret += atomName + "_";
        
        if (markerCount >= 1) {
            ret += to_string(edgeCount(atomIndex));
        } else {
            ret += "*";
        }
        
        if (markerCount >= 2) {
            ret += markAtomLink(atomIndex);
        } else {
            ret += "*";
        }
        
        if (markerCount >= 3) {
            ret += atomType(atomIndex);
        } else {
            ret += "*";
        }
        
        return ret;
    }        
    
    string chainName(const vector<int>& chain) {
        string ret;
        for (auto it : chain) {
            ret += markAtom(it);
        }
        if (ret.substr(0, 5) > ret.substr(ret.size() - 5, 5)) {
            ret = ret.substr(ret.size() - 5, 5) + ret.substr(5, ret.size() - 10) + ret.substr(0, 5);
        }
        return ret;
    }
    
    static vector<string> uniqueList(vector<string> names) {
        vector<vector<string>> decomposed(names.size());
        for (unsigned i = 0; i < names.size(); i++) {
            for (unsigned k = 0; k < names[i].size() / 5; k++) {
                decomposed[i].push_back(names[i].substr(k * 5, 5));
            }
        }
        for (int i = 0; i < int(decomposed.size()); i++) {
            for (int k = 0; k < i; k++) {
                if (equal(decomposed[i].begin(), decomposed[i].end(), decomposed[k].rbegin())) {
                    decomposed[k] = decomposed[i];
                }
            }
        }
        vector<string> ret(decomposed.size());
        for (unsigned i = 0; i < decomposed.size(); i++) {
            ret[i] = accumulate(decomposed[i].begin(), decomposed[i].end(), string(""));
        }
        return ret;
    }
    
    // создаёт словарь цепочек
    map<string, int> createList(int K) {
        vector<string> names;
        Graph g = createGraph();
        vector<vector<int>> chains = g.allSubgraphs(K);
        for (auto& it : chains) {
            names.push_back(chainName(it));
        }
        
        names = uniqueList(names);
        
        map<string, int> ret;
        for (auto& name : names) {
            if (ret.count(name) == 0) {
                ret[name] = 1;
            } else {
                ret[name]++;
            }
        }
        return ret;
    }
    
    char atomType(int vertex) {
        Graph g = createGraph();
        vector<int> reachable = g.reachableVertices(vertex);  
        int matchesCount = 0;
        int notMatchesCount = 0;
        for (int e = 0; e < g.vertexCount; e++) {
            if (e == vertex) continue;
            Graph copy = g;
            if (g.table[vertex][e] != 0 || g.table[e][vertex] != 0) {
                copy.table[vertex][e] = copy.table[e][vertex] = 0;
                if (copy.reachableVertices(vertex) != reachable) {
                    notMatchesCount++;
                } else {
                    matchesCount++;
                }
            }
        }
        if (matchesCount == 0) {
            return 'c';
        } else if (notMatchesCount == 0) {
            return 'r';
        } else {
            return 's';
        }            
    }
};

ostream& operator<<(ostream& out, Molecule mol) {
    out << "name=" << mol.name << ", atomCount=" << mol.atomCount << endl;
    for (auto it : mol.atoms) {
        out << "\t" << it << endl;
    }
    for (auto it : mol.links) {
        out << "\t" << it << endl;
    }
    return out;
}

vector<Molecule> loadStr(string filename) {
    vector<Molecule> ret;
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "loadStr: can\'t open the file \"" << filename << "\"";
        exit(-1);
    }
    string line;
    while (getline(file, line)) {
        stringstream ss(line);
        Molecule mol;
        ss >> mol.name >> mol.atomCount;
        for (int i = 0; i < mol.atomCount; i++) {
            getline(file, line);
            Atom a;
            int _;
            stringstream ss(line);
            ss >> a.index >> a.name;
            if (a.name.size() > 2) {
                a.hydrogenCount = stoi(a.name.substr(2));
                a.name = a.name.substr(0, 2);
            } else {
                ss >> a.hydrogenCount;
            }
            a.index--;
            ss >> _;
            for (int i = 0; i < 6; i++) {
                string linkAsStr;
                ss >> linkAsStr;
                if (stoi(linkAsStr) != 0) {
                    Link l;
                    l.fst = a.index;
                    l.snd = stoi(linkAsStr.substr(0, linkAsStr.size() - 1)) - 1;
                    if (l.fst < l.snd) {    // чтобы избежать повторов
                        l.type = stoi(linkAsStr.substr(linkAsStr.size() - 1));
                        if (l.type >= 5) {
                            l.type -= 4;
                        }
                        mol.links.push_back(l);
                    }
                }
            }
            ss >> a.x >> a.y >> a.z;
            mol.atoms.push_back(a);
        }
        ret.push_back(mol);
    }
    return ret;
}

vector<Molecule> loadSdf(string filename) {
    vector<Molecule> ret;
    ifstream file(filename);
    if (!file.is_open()) {
        cout << "loadSdf: can\'t open the file \"" << filename << "\"";
        exit(-1);
    }
    string line;
    stringstream ss;
    while (!file.eof()) {
        Molecule mol;
        int bondCount;
        getline(file, line);
        if (line == "") {
            break;
        }
        mol.name = line;
        getline(file, line);
        getline(file, line);
        getline(file, line);
        ss.str(line);
        ss >> mol.atomCount >> bondCount;
        for (int i = 0; i < mol.atomCount; i++) {
            Atom a;
            a.hydrogenCount = 0;
            getline(file, line);
            ss.str(line);
            ss >> a.x >> a.y >> a.z >> a.name;
            a.index = i;
            mol.atoms.push_back(a);
        }
        for (int i = 0; i < bondCount; i++) {
            Link l;
            getline(file, line);
            ss.str(line);
            ss >> l.fst >> l.snd >> l.type;
            l.fst--;
            l.snd--;
            mol.links.push_back(l);
        }
        ret.push_back(mol);
        while (getline(file, line) && line != "$$$$");
    }
    return ret;
}

map<string, int> createAllChains(vector<Molecule> mols, int K) {
    map<string, int> ret;
    for (auto it : mols) {
        auto list = it.createList(K);
        for (auto it : list) {
            if (ret.find(it.first) != ret.end()) {
                ret[it.first] += it.second;
            } else {
                ret[it.first] = it.second;
            }
        }
    }
    return ret;
}

//                      цепочка
// молекула из файла    *кол-во этих цепочек*
vector<int> createTable(vector<Molecule> mols, int K) {
    auto chains = createAllChains(mols, K);
    vector<int> ret(mols.size() * chains.size(), 0);
    for (unsigned i = 0; i < mols.size(); i++) {
        auto list = mols[i].createList(K);
        int chainIndex = 0;
        for (auto chain : chains) {
            for (auto it : list) {
                if (it.first == chain.first) {
                    ret[i * chains.size() + chainIndex] = it.second;
                }
            }
            chainIndex++;
        }
    }
    return ret;
}

void loadTableToFile(vector<Molecule> mols, int K, ostream& file) {
    auto table = createTable(mols, K);
    auto chains = createAllChains(mols, K);
    int i = 0;
    for (auto it : table) {
        file << it << " ";
        if (i == int(chains.size() - 1)) {
            file << endl;
            i = 0;
        } else {
            i++;
        }
    }
}

vector<Molecule> load(string filename) {
    string format = filename.substr(filename.rfind('.') + 1);
    if (format == "sdf" || format == "mol") {
        return loadSdf(filename);
    } else if (format == "str") {
        return loadStr(filename);
    } else {
        cout << "load: unknown file format \"" << filename << "\"";
        exit(-1);
    }
    return {};
}

template <typename Key, typename T>
void writeMapVert(std::ostream& out, const map<Key, T>& map) {
    for (auto& it : map) {
        out << it.first << ": " << it.second << "\n";
    }
}

class MolFiles {
public:
    string filename;
    vector<Molecule> mols;
    map<string, int> allChainsk2m1;
    map<string, int> allChainsk3m1;
    map<string, int> allChainsk2m2;
    map<string, int> allChainsk3m2;
    map<string, int> allChainsk2m3;
    map<string, int> allChainsk3m3;

    MolFiles(string _filename) {
        filename = "folder";
        mols = load(_filename);
		
        markerCount = 1; allChainsk2m1 = createAllChains(mols, 2);
        markerCount = 1; allChainsk3m1 = createAllChains(mols, 3);
        markerCount = 2; allChainsk2m2 = createAllChains(mols, 2);
        markerCount = 2; allChainsk3m2 = createAllChains(mols, 3);
        markerCount = 3; allChainsk2m3 = createAllChains(mols, 2);
        markerCount = 3; allChainsk3m3 = createAllChains(mols, 3);
    }
	
	void saveAllChains() {
		ofstream _("_");
        ofstream chainsk2m1File(filename + "/allChainsk2m1.txt");
        ofstream chainsk3m1File(filename + "/allChainsk3m1.txt");
        ofstream chainsk2m2File(filename + "/allChainsk2m2.txt");
        ofstream chainsk3m2File(filename + "/allChainsk3m2.txt");
        ofstream chainsk2m3File(filename + "/allChainsk2m3.txt");
        ofstream chainsk3m3File(filename + "/allChainsk3m3.txt");
		
		LOG(chainsk2m1File.is_open());
		LOG(chainsk3m1File.is_open());
		LOG(chainsk2m2File.is_open());
		LOG(chainsk3m2File.is_open());
		LOG(chainsk2m3File.is_open());
		LOG(chainsk3m3File.is_open());
		
//		LOG(mols);
		
  //      try{
		chainsk2m1File << allChainsk2m1;
		chainsk3m1File << allChainsk3m1;
		chainsk2m2File << allChainsk2m2;
		chainsk3m2File << allChainsk3m2;
		chainsk2m3File << allChainsk2m3;
		chainsk3m3File << allChainsk3m3;
		
		chainsk2m1File.close();
		chainsk3m1File.close();
		chainsk2m2File.close();
		chainsk3m2File.close();
		chainsk2m3File.close();
		chainsk3m3File.close();

    }
    
    void saveMolChains() {
		ofstream _("_");
        ofstream chainsk2m1File(filename + "/molVertChainsk2m1.txt");
        ofstream chainsk3m1File(filename + "/molVertChainsk3m1.txt");
        ofstream chainsk2m2File(filename + "/molVertChainsk2m2.txt");
        ofstream chainsk3m2File(filename + "/molVertChainsk3m2.txt");
        ofstream chainsk2m3File(filename + "/molVertChainsk2m3.txt");
        ofstream chainsk3m3File(filename + "/molVertChainsk3m3.txt");
		
		LOG(chainsk2m1File.is_open());
		LOG(chainsk3m1File.is_open());
		LOG(chainsk2m2File.is_open());
		LOG(chainsk3m2File.is_open());
		LOG(chainsk2m3File.is_open());
		LOG(chainsk3m3File.is_open());
        
		//LOG(mols);
		LOG(mols.size());
        
		for (Molecule mol : mols) {
        	chainsk2m1File << mol.name << ":\n";
	        chainsk3m1File << mol.name << ":\n";
	        chainsk2m2File << mol.name << ":\n";
	        chainsk3m2File << mol.name << ":\n";
	        chainsk2m3File << mol.name << ":\n";
	        chainsk3m3File << mol.name << ":\n";
	        
	        markerCount = 1; writeMapVert(chainsk2m1File, mol.createList(2));
	        markerCount = 1; writeMapVert(chainsk3m1File, mol.createList(3));
	        markerCount = 2; writeMapVert(chainsk2m2File, mol.createList(2));
	        markerCount = 2; writeMapVert(chainsk3m2File, mol.createList(3));
	        markerCount = 3; writeMapVert(chainsk2m3File, mol.createList(2));
	        markerCount = 3; writeMapVert(chainsk3m3File, mol.createList(3));
        }
    }
    
    void saveMatrices() {
		ofstream _("_");
        ofstream matrk2m1File(filename + "/matrk2m1.txt");
        ofstream matrk3m1File(filename + "/matrk3m1.txt");
        ofstream matrk2m2File(filename + "/matrk2m2.txt");
        ofstream matrk3m2File(filename + "/matrk3m2.txt");
        ofstream matrk2m3File(filename + "/matrk2m3.txt");
        ofstream matrk3m3File(filename + "/matrk3m3.txt");
		
		LOG(matrk2m1File.is_open());
		LOG(matrk3m1File.is_open());
		LOG(matrk2m2File.is_open());
		LOG(matrk3m2File.is_open());
		LOG(matrk2m3File.is_open());
		LOG(matrk3m3File.is_open());
        
        markerCount = 1; loadTableToFile(mols, 2, matrk2m1File); matrk2m1File << "\n";
        markerCount = 1; loadTableToFile(mols, 3, matrk3m1File); matrk3m1File << "\n";
        markerCount = 2; loadTableToFile(mols, 2, matrk2m2File); matrk2m2File << "\n";
        markerCount = 2; loadTableToFile(mols, 3, matrk3m2File); matrk3m2File << "\n";
        markerCount = 3; loadTableToFile(mols, 2, matrk2m3File); matrk2m3File << "\n";
        markerCount = 3; loadTableToFile(mols, 3, matrk3m3File); matrk3m3File << "\n";
    }
    
    void save(){
        system((string("mkdir ") + filename).c_str());
        saveAllChains();
        //saveMolChains();
        //saveMatrices();
    }
};

int main() {
    MolFiles molFiles("GLASS.sdf");
    molFiles.save();
    
    return 0;
}
