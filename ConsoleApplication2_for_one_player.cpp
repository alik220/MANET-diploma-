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

	const int player1AgentsNumber = 10;

	/*
		It's much more convenient to fill up
		a 2-D array than to call constructor
		functions a lot of times with different arguments
	*/
	int player1Pos[player1AgentsNumber][2] = {
		{1,1},{2,1},{2,2},
		{2,3},{2,4},{3,4},
		{4,4},{4,3},{4,2},
		{5,2}
	};

	Agent** player1Agents = new Agent * [player1AgentsNumber];

	for (int i = 0; i < player1AgentsNumber; i++)
	{
		player1Agents[i] = new Agent(player1Pos[i][0], player1Pos[i][1], 0x01);
		network.Add(player1Agents[i]->position, player1Agents[i]);
	}

	
	int player1Diameter = network.GetDiameter(
		[](Agent x) {return (x.playerID == 0x01); },
		[](Agent x) {return (x.playerID == 0x01); },
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


	vector<Pos> candidates;

	for (set<Pos, PosComp>::iterator iter = candidatesUniquifier.begin(); iter != candidatesUniquifier.end(); iter++)
	{
		candidates.push_back(*iter);
	}

	CalculateDistance distanceCalculator;

	Agent* player1Drone = new Agent(0x11);

	map<int, Pos> paretoOptimal;

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

			bool connected = (idist + jdist) == 1;
			bool haveNeighbors = player1NeighborsNumber > 0;
			bool haveMoreThanOneNeighbor = player1NeighborsNumber > 1;

			if ((haveNeighbors && connected) || haveMoreThanOneNeighbor)
			{
				player1Drone->position = candidates[i];
				network.Add(candidates[i], player1Drone);

				int player1Gain =
					player1Diameter - network.GetDiameter(
						[](Agent x) {return ((x.playerID == 0x01) || ((x.playerID & 0x10) != 0)); },
						[](Agent x) {return (x.playerID == 0x01); },
						(INT_MAX - 1) / 2
					);


				/*
					At this point we have the drones set
					and diameter calculated so we can perform
					any computations needed to figure out if
					this set of positions is valuable or not
				*/

				// seeking Pareto optimal solutions
				//pair<int, int> newSolution = make_pair(player1Gain, player2Gain);

				bool goesInside = true;

				for (
					map<int, Pos>::iterator iter = paretoOptimal.begin();
					iter != paretoOptimal.end();
					)
				{
					int checkAgainstPlayer1Gain = iter->first;
				}

			
				network.Remove(candidates[i], player1Drone);
			}
		}
	}

	double bestPageRank = 0.0;
	pair<Pos, Pos> solution;

	auto filter = [](Agent) {return true; };
	for (map<int, Pos>::iterator iter = paretoOptimal.begin(); iter != paretoOptimal.end(); iter++)
	{
		for (Pos::iterator setiter = iter->begin(); setiter != iter->second.end(); setiter++)
		{
			Pos p1pos = paretoOptimal;
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
