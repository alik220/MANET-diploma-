// ConsoleApplication2.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include "Header.h"

#include <map>
#include <set>
#include <string>
#include <iostream>
#include <vector>
#include <climits>

using namespace std;

class Pos
{
public:
	int i;
	int j;

	Pos() : i(0), j(0) {}

	Pos(int iInitial, int jInitial) : i(iInitial), j(jInitial) { }

	Pos& operator=(Pos rhs)
	{
		i = rhs.i;
		j = rhs.j;
		return *this;
	}

	Pos(const Pos& another) : i(another.i), j(another.j) {}

	Pos Left() { return Pos(this->i - 1, this->j); }
	Pos Right() { return Pos(this->i + 1, this->j); }
	Pos Up() { return Pos(this->i, this->j - 1); }
	Pos Down() { return Pos(this->i, this->j + 1); }

	bool operator<(const Pos another) const
	{
		return(this->i < another.i) || ((this->i == another.i) && (this->j < another.j));
	}
};

class PosComp
{
public:
	bool operator()(const Pos& lhs, const Pos& rhs) const
	{
		return (lhs.i < rhs.i) || ((lhs.i == rhs.i) && (lhs.j < rhs.j));
	}
};

class Agent
{
public:
	Pos position;

	int playerID;

	~Agent() {}

	Agent() : position(), playerID(-1) {}

	Agent(int iInitial, int jInitial, int playerID) : position(iInitial, jInitial), playerID(playerID) {}

	Agent(Pos positionInitial, int playerID) : position(positionInitial), playerID(playerID) {}

	Agent(int playerID) : playerID(playerID) {}

	Agent(Agent& another) : position(another.position), playerID(another.playerID) {}
};

class DoNeighbor
{
public:
	bool operator()(const Agent& lhs, const Agent& rhs)
	{
		int idist = abs(lhs.position.i - rhs.position.i);
		int jdist = abs(lhs.position.j - rhs.position.j);

		bool distanceIsRight = ((idist + jdist) <= 1) && ((idist + jdist) > 0);
		bool eitherIsADrone = ((lhs.playerID | rhs.playerID) & 0x10) != 0;
		bool samePlayer = (lhs.playerID == rhs.playerID);

		//return (distanceIsRight && (eitherIsADrone || samePlayer));
		return distanceIsRight;
	}
};

class CalculateDistance
{
public:
	int operator()(const Agent& lhs, const Agent& rhs)
	{
		int idist = abs(lhs.position.i - rhs.position.i);
		int jdist = abs(lhs.position.j - rhs.position.j);
		return (idist + jdist);
	}
};

template<class _T>
class IBase
{
public:
	virtual ~IBase() {};
	virtual _T GetData() = 0;
	virtual bool Comp(IBase& another) = 0;
};

template<class _T>
class Inherited : public IBase<_T>
{
public:
	~Inherited() {}

	_T GetData()
	{
		return 4;
	}

	bool Comp(IBase<_T>& another)
	{
		_T data = this->GetData();
		_T anotherdata = another.GetData();
		return data == anotherdata;
	}
};



int main()
{
	Network<Pos, Agent, DoNeighbor, int, CalculateDistance, PosComp> network;

	const int player1AgentsNumber = 9;
	const int player2AgentsNumber = 5;


	/*
		It's much more convenient to fill up
		a 2-D array than to call constructor
		functions a lot of times with different arguments
	*/

	
	//координаты для озера с Владом 2 игрок
	int player1Pos[player1AgentsNumber][2] =
	{
		{2, 1}, {2, 2}, {1, 2}, {2, 3}, {1, 3},
		{1, 4}, {1, 5}, {2,5}, {3, 5}
	};

	
	/*int player2Pos[player2AgentsNumber][2] =
	{
		{4, 3}, {5, 3}, {5, 2}, {6, 2}, {7, 2},
		{7, 3}, {7, 4}, {7, 5}, {6, 5}, {5, 5}
	};*/
	
	//
	int player2Pos[player2AgentsNumber][2] =
	{
		{6, 5}, {7, 5}, {7, 4}, {7, 3}, {6, 3}
	}; 

	/*
	int player1Pos[player1AgentsNumber][2] = {
	{1,1},{2,1},{2,2},
	{2,3},{2,4},{3,4},
	{4,4},{4,3},{4,2},
	{5,2}
	};

	int player2Pos[player2AgentsNumber][2] = {
		{3,5},{4,5},{5,5},
		{6,5},{4,4},{4,3},
		{4,2},{5,2},{6,2},
		{7,2},{8,2}, {2,5},
		{1,5}
	};
	*/


/*	int player1Pos[player1AgentsNumber][2] = {
		{1,4}, {2,1}, {2,2}, {2,3}, {2,4},
		{2, 5}, {2,6}, {2,7}, {3,1}, {3,7},
		{4,1}, {4,2}, {4,3}, {4,5}, {4,6}, {4,7},
		{5,3}, {5,5}
	};
*/


	/*
	int player1Pos[player1AgentsNumber][2] =
	{
		{1, 1}, {1, 2}, {1, 3}, {2, 1}, {2, 2}
	};

	int player2Pos[player2AgentsNumber][2] =
	{
		{1, 1}, {1, 2}, {1, 3}, {2, 1}, {2, 2}
	};
	*/

	//int* player1PosPtr = &player1Pos[0][0];  // указатель на позиции первого игрока






	/*int player1Pos[player1AgentsNumber][2] = {
		{9,1},{9,2},{9,3},
		{9,4},{9,5},{10,4},
		{11,4},{12,4},{13,4},
		{14,4},{14,5},{14,6},
		{14,7},{14,8},{13,8},
		{12,8},{12,9},{12,10},
		{12,11},{11,11},{10,11},
		{9,11},{9,10},{9,9},
		{9,8},{8,11},{9,12},
		{7,11},{6,11},//{6,10},
		{7,10},{7,9},//{6,9},
		//{7,7}, {7,6},
		{6,7},
		{6,6},{6,5},{6,4},
		{8,3},{8,4},{8,5},
		{7,4}
	};*/

	

	/*int player2Pos[player2AgentsNumber][2] = {
		{9,1},{9,2},{9,3},
		{9,4},{9,5},{10,4},
		{11,4},{12,4},{13,4},
		{14,4},{14,5},{14,6},
		{14,7},{14,8},{13,8},
		{12,8},{12,9},{12,10},
		{12,11},{11,11},{10,11},
		{9,11},{9,10},{9,9},
		{9,8},{8,11},{8,12},
		{7,11},{6,11},{6,10},
		{7,10},{7,9},{6,9},
		{7,7},{7,6},{6,7},
		{6,6},{6,5},{6,4},
		{8,3},{8,4},{8,5},
		{7,4}
	};*/





	/*int player2Pos[player2AgentsNumber][2];

	for (int i = 0; i < player2AgentsNumber; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			player2Pos[i][j] = player1Pos[i][j];
		}
	}
	*/
	

	Agent** player1Agents = new Agent * [player1AgentsNumber];
	Agent** player2Agents = new Agent * [player2AgentsNumber];

	for (int i = 0; i < player1AgentsNumber; i++)
	{
		player1Agents[i] = new Agent(player1Pos[i][0], player1Pos[i][1], 0x01);
		network.Add(player1Agents[i]->position, player1Agents[i]);
	}

	for (int i = 0; i < player2AgentsNumber; i++)
	{
		player2Agents[i] = new Agent(player2Pos[i][0], player2Pos[i][1], 0x02);
		network.Add(player2Agents[i]->position, player2Agents[i]);
	}

	int player1Diameter = network.GetDiameter(
		[](Agent x) {return (x.playerID == 0x01); },
		[](Agent x) {return (x.playerID == 0x01); },
		(INT_MAX - 1) / 2
	);

	int player2Diameter = network.GetDiameter(
		[](Agent x) {return (x.playerID == 0x02); },
		[](Agent x) {return (x.playerID == 0x02); },
		(INT_MAX - 1) / 2
	);

	set<Pos, PosComp> existingPositions;
	set<Pos, PosComp> candidatesUniquifier;

	for (int i = 0; i < player1AgentsNumber; i++)
	{
		candidatesUniquifier.insert(Pos(player1Pos[i][0] - 1, player1Pos[i][1]));
		candidatesUniquifier.insert(Pos(player1Pos[i][0] + 1, player1Pos[i][1]));
		candidatesUniquifier.insert(Pos(player1Pos[i][0], player1Pos[i][1] - 1));
		candidatesUniquifier.insert(Pos(player1Pos[i][0], player1Pos[i][1] + 1));

		existingPositions.insert(Pos(player1Pos[i][0], player1Pos[i][1]));
	}

	for (int i = 0; i < player2AgentsNumber; i++)
	{
		candidatesUniquifier.insert(Pos(player2Pos[i][0] - 1, player2Pos[i][1]));
		candidatesUniquifier.insert(Pos(player2Pos[i][0] + 1, player2Pos[i][1]));
		candidatesUniquifier.insert(Pos(player2Pos[i][0], player2Pos[i][1] - 1));
		candidatesUniquifier.insert(Pos(player2Pos[i][0], player2Pos[i][1] + 1));

		existingPositions.insert(Pos(player2Pos[i][0], player2Pos[i][1]));
	}

	vector<Pos> candidates;

	for (set<Pos, PosComp>::iterator iter = candidatesUniquifier.begin(); iter != candidatesUniquifier.end(); iter++)
	{
		candidates.push_back(*iter);
	}

	CalculateDistance distanceCalculator;

	Agent* player1Drone = new Agent(0x11);
	Agent* player2Drone = new Agent(0x12);

	map<pair<int, int>, set<pair<Pos, Pos>>> paretoOptimal;

	for (int i = 0; i < candidates.size(); i++)
	{
		for (int j = i; j < candidates.size(); j++)
		{
			/*
				checking if the pair of positions
				can possibly provide some improvement
			*/

			int idist = abs(candidates[i].i - candidates[j].i);
			int jdist = abs(candidates[i].j - candidates[j].j);

			int player1NeighborsNumber =
				(existingPositions.find(candidates[i].Left()) != existingPositions.end()) +
				(existingPositions.find(candidates[i].Right()) != existingPositions.end()) +
				(existingPositions.find(candidates[i].Up()) != existingPositions.end()) +
				(existingPositions.find(candidates[i].Down()) != existingPositions.end());

			int player2NeighborsNumber =
				(existingPositions.find(candidates[j].Left()) != existingPositions.end()) +
				(existingPositions.find(candidates[j].Right()) != existingPositions.end()) +
				(existingPositions.find(candidates[j].Up()) != existingPositions.end()) +
				(existingPositions.find(candidates[j].Down()) != existingPositions.end());

			bool connected = (idist + jdist) == 1;
			bool haveNeighbors = ((player1NeighborsNumber > 0) && (player2NeighborsNumber > 0));
			bool haveMoreThanOneNeighbor = ((player1NeighborsNumber > 1) && (player2NeighborsNumber > 1));

			if ((haveNeighbors && connected) || haveMoreThanOneNeighbor)
			{
				player1Drone->position = candidates[i];
				player2Drone->position = candidates[j];
				network.Add(candidates[i], player1Drone);
				network.Add(candidates[j], player2Drone);

				int player1Gain =
					player1Diameter - network.GetDiameter(
						[](Agent x) {return ((x.playerID == 0x01) || ((x.playerID & 0x10) != 0)); },
						[](Agent x) {return (x.playerID == 0x01); },
						(INT_MAX - 1) / 2
					);

				int player2Gain =
					player2Diameter - network.GetDiameter(
						[](Agent x) {return ((x.playerID == 0x02) || ((x.playerID & 0x10) != 0)); },
						[](Agent x) {return (x.playerID == 0x02); },
						(INT_MAX - 1) / 2
					);

				/*
					At this point we have the drones set
					and diameter calculated so we can perform
					any computations needed to figure out if
					this set of positions is valuable or not
				*/

				// seeking Pareto optimal solutions
				pair<int, int> newSolution = make_pair(player1Gain, player2Gain);

				bool goesInside = true;

				for (
					map<pair<int, int>, set<pair<Pos, Pos>>>::iterator iter = paretoOptimal.begin();
					iter != paretoOptimal.end();
					)
				{
					int checkAgainstPlayer1Gain = iter->first.first;
					int checkAgainstPlayer2Gain = iter->first.second;

					if (player1Gain == checkAgainstPlayer1Gain && player2Gain == checkAgainstPlayer2Gain)
					{
						iter++;
						continue;
					}

					if (player1Gain <= checkAgainstPlayer1Gain && player2Gain <= checkAgainstPlayer2Gain)
					{
						goesInside = false;
						break;
					}

					if (!(player1Gain >= checkAgainstPlayer1Gain && player2Gain >= checkAgainstPlayer2Gain))
					{
						iter++;
					}
					else
					{
						iter = paretoOptimal.erase(iter);
					}
				}

				if (goesInside)
				{
					paretoOptimal[newSolution].insert(make_pair(candidates[i], candidates[j]));
				}

				network.Remove(candidates[i], player1Drone);
				network.Remove(candidates[j], player2Drone);
			}
		}
	}
	

	
	double bestBetwCentr = 0.0;
	pair<Pos, Pos> solution2;
	vector<pair<Agent*, double>> betwCentr;

	auto filter1 =  [](Agent x) {
		if (x.playerID == 0x01)
			return true;
		else
			return false;
 } ;
	auto filter2 = [](Agent x) {
		if (x.playerID == 0x02)
			return true;
		else
			return false;
	};
	//auto filter2 = [](Agent) {return true; };
	vector< vector <int> > temp;
	temp = network.GreedyAl(filter1);
	for (map<pair<int, int>, set<pair<Pos, Pos>>>::iterator iter = paretoOptimal.begin(); iter != paretoOptimal.end(); iter++)
	{
		for (set<pair<Pos, Pos>>::iterator setiter = iter->second.begin(); setiter != iter->second.end(); setiter++)
		{
			Pos p1pos = setiter->first;
			Pos p2pos = setiter->second;

			player1Drone->position = p1pos;
			player2Drone->position = p2pos;
			network.Add(p1pos, player1Drone);
			network.Add(p2pos, player2Drone);

			betwCentr = network.BetwCentral(filter1);

			double currentBetwCentr = 0.0;

			for (int i = 0; i < betwCentr.size(); i++)
			{
				if (betwCentr[i].first == player1Drone || betwCentr[i].first == player2Drone)
				{
					currentBetwCentr += betwCentr[i].second;
					///cout << currentBetwCentr << endl;
				}
			}

			if (currentBetwCentr > bestBetwCentr)
			{
				bestBetwCentr = currentBetwCentr;
				solution2 = make_pair(p1pos, p2pos);
			}

			network.Remove(p1pos, player1Drone);
			network.Remove(p2pos, player2Drone);
		}
	}

	delete player1Drone;
	delete player2Drone;

	for (int i = 0; i < player1AgentsNumber; i++)
	{
		network.Remove(player1Agents[i]->position, player1Agents[i]);
		delete player1Agents[i];
	}

	for (int i = 0; i < player2AgentsNumber; i++)
	{
		network.Remove(player2Agents[i]->position, player2Agents[i]);
		delete player2Agents[i];
	}

	delete[] player1Agents;
	delete[] player2Agents;

	cout << "(" << solution2.first.i << ", " << solution2.first.j << ") (" <<
		solution2.second.i << ", " << solution2.second.j << ")" << "; "<< "Best betweenes centrality: " << bestBetwCentr << endl;
	




	/*
	double bestPageRank = 0.0;
	pair<Pos, Pos> solution;

	auto filter = [](Agent) {return true; };
	for (map<pair<int, int>, set<pair<Pos, Pos>>>::iterator iter = paretoOptimal.begin(); iter != paretoOptimal.end(); iter++)
	{
		for (set<pair<Pos, Pos>>::iterator setiter = iter->second.begin(); setiter != iter->second.end(); setiter++)
		{
			Pos p1pos = setiter->first;
			Pos p2pos = setiter->second;

			player1Drone->position = p1pos;
			player2Drone->position = p2pos;
			network.Add(p1pos, player1Drone);
			network.Add(p2pos, player2Drone);

			vector<pair<Agent*, double>> pageRank = network.GetPageRank(filter, filter);

			double currentPageRank = 0.0;

			for (int i = 0; i < pageRank.size(); i++)
			{
				if (pageRank[i].first == player1Drone || pageRank[i].first == player2Drone)
				{
					currentPageRank += pageRank[i].second;
				}
			}

			if (currentPageRank > bestPageRank)
			{
				bestPageRank = currentPageRank;
				solution = make_pair(p1pos, p2pos);
			}

			network.Remove(p1pos, player1Drone);
			network.Remove(p2pos, player2Drone);
		}
	}

	delete player1Drone;
	delete player2Drone;

	for (int i = 0; i < player1AgentsNumber; i++)
	{
		network.Remove(player1Agents[i]->position, player1Agents[i]);
		delete player1Agents[i];
	}

	for (int i = 0; i < player2AgentsNumber; i++)
	{
		network.Remove(player2Agents[i]->position, player2Agents[i]);
		delete player2Agents[i];
	}

	delete[] player1Agents;
	delete[] player2Agents;

	cout << "(" << solution.first.i << ", " << solution.first.j << ") (" <<
		solution.second.i << ", " << solution.second.j << ")" << endl;

	*/






//расчет всех оптимальных решений без мер центральности

/*
pair<Pos, Pos> solution;

auto filter = [](Agent) {return true; };
	for (map<pair<int, int>, set<pair<Pos, Pos>>>::iterator iter = paretoOptimal.begin(); iter != paretoOptimal.end(); iter++)
	{
		for (set<pair<Pos, Pos>>::iterator setiter = iter->second.begin(); setiter != iter->second.end(); setiter++)
		{
			Pos p1pos = setiter->first;
			Pos p2pos = setiter->second;

			player1Drone->position = p1pos;
			player2Drone->position = p2pos;
			network.Add(p1pos, player1Drone);
			network.Add(p2pos, player2Drone);

			solution = make_pair(p1pos, p2pos);
			
			//player1Diameter = betw_cent(solution.first.i, solution.first.j,
			//solution.second.i, solution.second.j, player1PosPtr, player1AgentsNumber);
			
			cout << "(" << solution.first.i << ", " << solution.first.j << ")  (" <<
				solution.second.i << ", " << solution.second.j << ")" << endl;
			

		}
	}
	

	cout << "(" << solution.first.i << ", " << solution.first.j << ") (" <<
	solution.second.i << ", " << solution.second.j << ")" << endl; 
	
*/

	




	
	fgetc(stdin);

	return 0;
}


// Запуск программы: CTRL+F5 или меню "Отладка" > "Запуск без отладки"
// Отладка программы: F5 или меню "Отладка" > "Запустить отладку"

// Советы по началу работы 
//   1. В окне обозревателя решений можно добавлять файлы и управлять ими.
//   2. В окне Team Explorer можно подключиться к системе управления версиями.
//   3. В окне "Выходные данные" можно просматривать выходные данные сборки и другие сообщения.
//   4. В окне "Список ошибок" можно просматривать ошибки.
//   5. Последовательно выберите пункты меню "Проект" > "Добавить новый элемент", чтобы создать файлы кода, или "Проект" > "Добавить существующий элемент", чтобы добавить в проект существующие файлы кода.
//   6. Чтобы снова открыть этот проект позже, выберите пункты меню "Файл" > "Открыть" > "Проект" и выберите SLN-файл.
