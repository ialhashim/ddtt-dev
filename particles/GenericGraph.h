/* Retrieved from: http://en.literateprograms.org/Dijkstra's_algorithm_(C_Plus_Plus)?oldid=13422 */
#pragma once

#include <vector>
#include <map>
#include <list>
#include <queue>
#include <set>

#undef min
#undef max
#include <limits>

#define Max(a,b) (((a) > (b)) ? (a) : (b))
#define Min(a,b) (((a) < (b)) ? (a) : (b))

static size_t GLOBAL_GRAPH_UID = 0;

namespace GenericGraphs{

	template <typename VertexType = unsigned int, typename WeightType = double>
	class Graph
	{
	private:
		typedef VertexType vertex_t;
		typedef WeightType weight_t;

	public:
		struct Edge {
			vertex_t target;
			weight_t weight;
			unsigned int index;

			Edge(vertex_t arg_target, weight_t arg_weight, unsigned int edge_index = -1)
				: target(arg_target), weight(arg_weight), index(edge_index) { }

			bool operator == (const Edge & rhs) const
			{
				return (target == rhs.target && weight == rhs.weight && index == rhs.index);
			}
		};

		struct CompareEdge{
			bool operator()(const Edge &a, const Edge &b){
				if (a.target < b.target) return true;
				if (a.target > b.target) return false;
				return a.index < b.index;
			}
		};

		typedef std::set<Edge,CompareEdge> EdgesSet;

	private:
		typedef std::map<vertex_t, std::list<Edge> > adjacency_map_t;
		typedef std::set<vertex_t> vertices_set;

		template <typename T1, typename T2>
		struct pair_first_less{
			bool operator()(std::pair<T1,T2> p1, std::pair<T1,T2> p2) const {
				if(p1.first == p2.first) {
					return p1.second < p2.second;
				}
				return p1.first < p2.first;
			}
		};

	public:
		void DijkstraComputePathsMany(const std::set<vertex_t> & sources)
		{
			// Add dummy start node
			auto newNode = AddVertex( vertex_t(vertices.size()) );

			// Add edges to dummy node
			for(auto & v : sources)	AddEdge(v, newNode, weight_t(0), -1);

			DijkstraComputePaths( newNode );

			// Remove edges
			for(auto & v : sources) removeEdge(v, newNode);

			// Remove dummy node
			vertices.erase( newNode );
		}

		void DijkstraComputePaths(vertex_t source)
		{
			if(source == lastStart)
				return;
			
			lastStart = source;
			previous.clear();

			min_distance.clear();
			min_distance.resize(vertices.size(), std::numeric_limits< WeightType >::infinity());

			for (typename adjacency_map_t::iterator vertex_iter = adjacency_map.begin();
				vertex_iter != adjacency_map.end(); vertex_iter++){
					vertex_t v = vertex_iter->first;
					previous[v] = -1;
			}

			min_distance[source] = 0;

			std::set< std::pair<weight_t, vertex_t>, pair_first_less<weight_t, vertex_t> > vertex_queue;
			vertex_queue.insert(std::make_pair(min_distance[source], source));

			while (!vertex_queue.empty()) {
				vertex_t u = vertex_queue.begin()->second;
				vertex_queue.erase(vertex_queue.begin());

				// Visit each edge exiting u
				for (typename std::list<Edge>::iterator edge_iter = adjacency_map[u].begin();
					edge_iter != adjacency_map[u].end();
					edge_iter++)
				{
					vertex_t v = edge_iter->target;
					weight_t weight = edge_iter->weight;
					weight_t distance_through_u = min_distance[u] + weight;
					if (distance_through_u < min_distance[v]) {
						vertex_queue.erase(std::pair<weight_t, vertex_t>(min_distance[v], v));

						min_distance[v] = distance_through_u;
						previous[v] = u;
						vertex_queue.insert(std::pair<weight_t, vertex_t>(min_distance[v], v));
					}
				}
			}
		}

		std::list<vertex_t> DijkstraGetShortestPathsTo(vertex_t target)
		{
			std::list<vertex_t> path;
			typename std::map<vertex_t, vertex_t>::iterator prev;
			vertex_t vertex = target;

			path.push_front(vertex);

			while((prev = previous.find(vertex)) != previous.end())
			{
				vertex = prev->second;
				path.push_front(vertex);
			}

			// Remove last 'previous' which is null
			if(path.size() > 1) path.erase(path.begin());

			return path;
		}

		std::list<vertex_t> DijkstraShortestPath(vertex_t start, vertex_t end)
		{
			this->DijkstraComputePaths(start);
			return this->DijkstraGetShortestPathsTo(end);
		}

		std::vector<Edge> DijkstraShortestPathEdges(vertex_t start, vertex_t end)
		{
			this->DijkstraComputePaths(start);
			auto path = this->DijkstraGetShortestPathsTo(end);

			std::vector<Edge> result;
			if(path.size() < 2) return result;

			vertex_t lastVert = start;
			for(auto v : path){
				result.push_back(Edge(lastVert, 1.0, v));
				lastVert = v;
			}
			result.erase(result.begin());
			return result;
		}

		int NodeDistance(vertex_t n1, vertex_t n2)
		{
			if(n1 == n2) return 0;
			if(!this->isConnected(n1, n2)) return -1;
			return DijkstraShortestPath(n1, n2).size();
		}

	public:
		// Graph Variables:
		vertices_set vertices;
		adjacency_map_t adjacency_map;
		vertex_t lastStart;

	public:

		std::vector<weight_t> min_distance;
		std::map<vertex_t, vertex_t> previous;

		Graph()
		{
			lastStart = std::numeric_limits<vertex_t>::max();
			sid = 0;

			uid = GLOBAL_GRAPH_UID++;
		}

		Graph(const Graph& from)
		{
			this->adjacency_map = from.adjacency_map;
			this->vertices = from.vertices;
			this->min_distance = from.min_distance;
			this->previous = from.previous;
			this->lastStart = from.lastStart;
			this->sid = from.sid;

			uid = GLOBAL_GRAPH_UID++; 
		}

		void AddEdge(vertex_t p1, vertex_t p2, weight_t weight, int index = -1)
		{
			if( IsEdgeExists(p1, p2) ) return;

			adjacency_map[AddVertex(p1)].push_back(Edge(AddVertex(p2), weight, index));
			adjacency_map[AddVertex(p2)].push_back(Edge(AddVertex(p1), weight, index));
		}

		vertex_t AddVertex(vertex_t p)
		{
			vertices.insert(p);

			return p;
		}

		vertex_t FirstVertex()
		{
			if(vertices.empty()) return -1;
			return *vertices.begin();
		}

		void removeDirectedEdge(vertex_t p1, vertex_t p2)
		{
			std::list<Edge> * adj = &adjacency_map[p1];

			for(typename std::list<Edge>::iterator i = adj->begin(); i != adj->end(); i++)
			{
				Edge * e = &(*i);

				if(e->target == p2)
				{
					adj->remove(*e);
					return;
				}
			}
		}

		void removeEdge(vertex_t p1, vertex_t p2)
		{
			removeDirectedEdge(p1, p2);
			removeDirectedEdge(p2, p1);
		}

		void SetDirectedEdgeWeight(vertex_t p1, vertex_t p2, weight_t newWeight)
		{
			std::list<Edge> * adj = &adjacency_map[p1];

			for(typename std::list<Edge>::iterator i = adj->begin(); i != adj->end(); i++)
			{
				Edge * e = &(*i);

				if(e->target == p2)
				{
					e->weight = newWeight;
					return;
				}
			}
		}

		void SetEdgeWeight(vertex_t p1, vertex_t p2, weight_t newWeight)
		{
			SetDirectedEdgeWeight(p1, p2, newWeight);
			SetDirectedEdgeWeight(p2, p1, newWeight);
		}

		vertex_t GetRandomNeighbour(vertex_t p)
		{
			return adjacency_map[p].front().target;
		}

		std::vector<vertex_t> GetNeighbours(vertex_t p)
		{
			std::vector<vertex_t> neighbours;

			std::list<Edge> * adj = &adjacency_map[p];

			for(typename std::list<Edge>::iterator i = adj->begin(); i != adj->end(); i++)
			{
				Edge * e = &(*i);

				neighbours.push_back(e->target);
			}

			return neighbours;
		}

		vertex_t GetOtherNeighbour(vertex_t p, vertex_t q)
		{
			vertex_t n = p;

			std::list<Edge> * adj = &adjacency_map[p];
			for(typename std::list<Edge>::iterator i = adj->begin(); i != adj->end(); i++)
			{
				Edge * e = &(*i);

				if(e->target != q)
				{
					n = e->target;
					break;
				}
			}

			return n;
		}

		bool isCircular(vertex_t p)
		{
			std::set<vertex_t> visited;

			visited.insert(p);

			bool hasMore = true;

			vertex_t curr = p;
			vertex_t prev = p;

			while(hasMore)
			{
				vertex_t next = GetOtherNeighbour(curr, prev);

				if(next == curr)
					return false;

				if(visited.find(next) != visited.end())
					return true;

				prev = curr;
				curr = next;
			}

			return false;
		}

		bool isConnected(vertex_t v1, vertex_t v2)
		{
			this->DijkstraComputePaths(v1);

			if(min_distance[v2] != std::numeric_limits< WeightType >::infinity())
				return true;

			return false;
		}

		vertices_set GetNodes()
		{
			return vertices;
		}

		bool IsEmpty()
		{
			return vertices.empty();
		}

		bool IsHasVertex( vertex_t v )
		{
			return vertices.find(v) != vertices.end();
		}

		inline bool IsEdgeExists(vertex_t v1, vertex_t v2)
		{
			std::list<Edge> * adj = &adjacency_map[v1];

			for(typename std::list<Edge>::iterator i = adj->begin(); i != adj->end(); i++){
				Edge * e = &(*i);
				if(e->target == v2) return true;
			}

			return false;
		}

		// the 'index' of the edge will be replaced with index of a vertex
		std::vector<Edge> GetEdges()
		{
			std::vector<Edge> result;

			for(typename adjacency_map_t::iterator it = adjacency_map.begin(); it != adjacency_map.end(); it++)
			{
				vertex_t v1 = it->first;
				std::list<Edge> adj = it->second;

				for(typename std::list<Edge>::iterator i = adj.begin(); i != adj.end(); i++)
				{
					Edge e = *i;

					result.push_back(Edge(Min(e.target, v1), e.weight, Max(e.target, v1)));
				}
			}

			return result;
		}

		EdgesSet GetEdgesSet()
		{
			std::vector<Edge> allEdges = GetEdges();

			EdgesSet result;

			for(typename std::vector<Edge>::iterator it = allEdges.begin(); it != allEdges.end(); it++)
				result.insert(*it);

			return result;
		}

		size_t GetNumberEdges()
		{
			return GetEdgesSet().size();
		}

		std::vector<vertex_t> GetLeaves() const
		{
			std::vector<vertex_t> leaves;

			for(typename adjacency_map_t::const_iterator it = adjacency_map.begin(); it != adjacency_map.end(); it++)
			{
				if(it->second.size() < 2)
					leaves.push_back(it->first);
			}

			return leaves;
		}

		void explore(vertex_t seed, std::set<vertex_t> & explored)
		{
			std::queue<vertex_t> q;

			q.push(seed);
			explored.insert(seed);

			while(!q.empty())
			{
				vertex_t i = q.front();
				q.pop();

				std::list<Edge> * adj = &adjacency_map[i];
				for(typename std::list<Edge>::const_iterator it = adj->begin(); it != adj->end(); it++)
				{
					vertex_t j = it->target;

					// Check: is not visited ?
					if(explored.find(j) == explored.end()) 
					{
						explored.insert(j);
						q.push(j);
					}
				}
			}
		}

		std::vector<std::set<vertex_t> > GetConnectedComponents()
		{
			std::vector< std::set<vertex_t> > connectedComponents;
			std::set<vertex_t> unvisited;

			// fill unvisited set
			for(typename vertices_set::const_iterator it = vertices.begin(); it != vertices.end(); it++)
				unvisited.insert(*it);

			while(unvisited.size() > 1)
			{
				// Take first unvisited node
				vertex_t firstNode = *(unvisited.begin());

				// Explore its tree
				std::set<vertex_t> currVisits;
				currVisits.insert(firstNode);
				explore(firstNode, currVisits);

				// Add as a connected component
				connectedComponents.push_back(currVisits);

				// Remove from unvisited set
				for(typename std::set<vertex_t>::iterator it = currVisits.begin(); it != currVisits.end(); it++)
					unvisited.erase(*it);
			}

			return connectedComponents;
		}

		std::set<vertex_t> GetLargestConnectedComponent()
		{
			std::vector<std::set<vertex_t> > connectedComponents = GetConnectedComponents();

			// Find set with maximum number of nodes
			int maxConnectSize = -1, max_i = 0;
			for(int i = 0; i < (int)connectedComponents.size(); i++)
			{
				int currSize = connectedComponents[i].size();

				if(currSize > maxConnectSize){
					maxConnectSize = currSize;
					max_i = i;
				}
			}

			// Return maximum set
			return connectedComponents[max_i];
		}

		void subGraph(Graph & g, const std::set<vertex_t> & explored)
		{
			for(typename std::set<vertex_t>::const_iterator vi = explored.begin(); vi != explored.end(); vi++)
			{
				std::list<Edge> adj = g.adjacency_map[*vi];

				for(typename std::list<Edge>::iterator e = adj.begin(); e != adj.end(); e++)
					this->AddEdge(*vi, e->target, e->weight);
			}

			// Isolated nodes
			if(this->vertices.empty())
			{
				for(typename std::set<vertex_t>::const_iterator vi = explored.begin(); vi != explored.end(); vi++)
					this->AddVertex( *vi );
			}
		}

		std::vector< Graph <vertex_t,weight_t> > toConnectedParts()
		{
			std::vector< Graph <vertex_t,weight_t> > result;

			// Make a 'bitmap' of visited nodes
			std::map<vertex_t, bool> isVisited;
			for(typename std::set<vertex_t>::iterator it = vertices.begin(); it != vertices.end(); it++)
				isVisited[*it] = false;

			for(typename std::map<vertex_t,bool>::iterator i = isVisited.begin(); i != isVisited.end(); i++)
			{
				// Check if visited
				if(i->second)
					continue;

				vertex_t seed = i->first;

				std::set<vertex_t> explored;
				explore(seed, explored);

				// Add this new connected sub graph from exploration
				result.push_back(Graph<vertex_t, weight_t>());
				result.back().subGraph(*this, explored);

				// mark as visited the explored nodes
				for(typename std::set<vertex_t>::iterator vi = explored.begin(); vi != explored.end(); vi++)
					isVisited[*vi] = true;
			}

			return result;
		}

		std::list<vertex_t> GetLargestConnectedPath()
		{
			std::list<vertex_t> longestPath;

			std::set<vertex_t> seedSet = GetLargestConnectedComponent(); 
			std::vector<vertex_t> leaves = GetLeaves();

			if(!leaves.size())
			{
				leaves.push_back(*this->vertices.begin());
			}

			vertex_t seed = leaves.front();

			DijkstraComputePaths(seed);

			for(typename std::set<vertex_t>::iterator it = seedSet.begin(); it != seedSet.end(); it++)
			{
				std::list<vertex_t> curPath = DijkstraGetShortestPathsTo(*it);

				if(curPath.size() > longestPath.size())
					longestPath = curPath;
			}

			return longestPath;
		}

		bool CheckAdjacent(vertex_t v1, vertex_t v2)
		{
			if(v1 == v2) return true;

			for(typename std::list<Edge>::iterator e = adjacency_map[v1].begin(); e != adjacency_map[v1].end(); e++)
				if(e->target == v2) return true;

			return false;
		}

		std::map< std::string, std::set<int> > intSet;
		size_t uid, sid;

		// Warning: memory leak properties
		std::map<std::string, void*> property;
		template<typename T>
		void setProperty(std::string pname, const T & value){
			T * ptr = new T;
			(*ptr) = value;
			property[pname] = (void *) ptr;
		}
		template<typename T>
		T * getProperty(std::string pname){
			if(property.find(pname) == property.end()) return NULL;
			return (T*)property[pname];
		}
	};
}
