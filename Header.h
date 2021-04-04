#pragma once

#ifndef MANET_PROJECT_NETWORK_HEADER_DEFINED
#define MANET_PROJECT_NETWORK_HEADER_DEFINED

#include <map>
#include <set>
#include <functional>
#include <queue>
#include <stack>
#include <string>
#include <iostream>
#include <vector>
#include <climits>

using namespace std;

/*
	We imply that anything can be a position
	for an agent of a network and anything
	can be an agent. So both of them are templates -
	TPosition and TAgents respectively.

	Compare must define a total order on a set
	of all possible TPositions.

	DoNeighbor is supposed to be a class or a struct
	that overrides an invocation operator - that is "()".
	The invokation operator must take 2 TAgents and return
	bool - true if these agents can communicate in
	some way, and false otherwise.

	Measure represents a type of data used to calculate
	distances on the network. And CalculateDistance must
	be a struct or a class that overrides an invocation
	operator so that it takes 2 TAgents and returns a
	value of type Measure, indicating distance between
	the TAgents.
*/
template <
	class TPosition,
	class TAgent,
	typename DoNeighbor,
	class Measure,
	typename CalculateDistance,
	typename Compare = less< TPosition >
>
class Network
{
private:
	int agentsCount = 0;

	class NetworkAgentData
	{
	public:
		TAgent* agent;
		set<NetworkAgentData*> neighbors;
		set<NetworkAgentData*> reverseNeighbors;

		NetworkAgentData()
		{
			agent = NULL;
			neighbors = set<NetworkAgentData*>();
			reverseNeighbors = set<NetworkAgentData*>();
		}

		NetworkAgentData(TAgent* agentPtr) : agent(agentPtr)
		{
			neighbors = set<NetworkAgentData*>();
			reverseNeighbors = set<NetworkAgentData*>();
		}

		NetworkAgentData(NetworkAgentData& original)
		{
			this->agent = original.agent;
			this->neighbors = original.neighbors;
			this->reverseNeighbors = original.reverseNeighbors;
		}

		~NetworkAgentData()
		{
			neighbors.clear();
			reverseNeighbors.clear();
			agent = NULL;
		}
	};

	class NetworkAgentDataComparator
	{
	public:
		bool operator()(NetworkAgentData* lhs, NetworkAgentData* rhs) const
		{
			/*
				This is an absolutely awesome way to
				sort these classes inside STL containers,
				since we may have two agents with absolutely
				identical properties and location (position),
				but we may never have two identical agent
				pointers direct to different agents
			*/
			return ((int)lhs->agent) < ((int)rhs->agent);
		}
	};

	/*
	 * From the definition you can figure out that
	 * networkAgents maps positions (coordinates)
	 * to the network agents that are located in
	 * those positions. Note, that a single position
	 * may contain as many agents as needed.
	 *
	 * Compare is a functor-class that defines a
	 * linear order over a set of all possible positions.
	 *
	 * It's not std::set because
	 */
	map<TPosition, set<NetworkAgentData*, NetworkAgentDataComparator>, Compare> networkAgents;

	/*
	 * class Iterator - simplifies iterating through networkAgents
	 *
	 * If you'd like to iterate through the whole set of existing
	 * agents, which is stored in a map of position - set_of_agents
	 * pairs, you would need a nested loop where the outer loop
	 * iterator would be a map<...>::iterator and the inner loop
	 * iterator would be a set<...>::iterator. This class flattens
	 * those nested loops, making the iteration process possible
	 * with a single loop. The usage is very similar to stl iterators -
	 * you initialize a beginning iterator with a BeginIterator() method,
	 * then you advance it one step forward with a pre-increment
	 * or post-increment operator and then to figure out when to stop
	 * you compare the current iterator against the EndIterator().
	 * Accessing the element to which the iterator points is done
	 * by using either dereferencing operator: * or ->
	 *
	 * Defines:
	 *		operator== and operator!= - check if iterators
	 *			are equal or if they are not equal respectively.
	 *
	 *		operator++() and operator++(int) - pre-increment
	 *			and post-increment respectively. Move iterator
	 *			one step forward, returning either an old or
	 *			an updated value, correspondingly to the type
	 *			of increment used.
	 *
	 *		operator*()
	 */
	class Iterator
	{
	public:
		typename map<TPosition, set<NetworkAgentData*, NetworkAgentDataComparator>, Compare>::iterator mapIterator;
		typename set<NetworkAgentData*, NetworkAgentDataComparator>::iterator setIterator;

		typename map<TPosition, set<NetworkAgentData*, NetworkAgentDataComparator>, Compare>::iterator mapEnd;
		typename set<NetworkAgentData*, NetworkAgentDataComparator>::iterator currentSetEnd;

		bool operator==(const Iterator& rhs)
		{
			return((this->mapIterator == rhs.mapIterator) && ((this->mapIterator == this->mapEnd) || (this->setIterator == rhs.setIterator)));
		}

		bool operator!=(const Iterator& rhs)
		{
			return(!((*this) == rhs));
		}

		Iterator& operator++()
		{
			if (mapIterator == mapEnd)
			{
				throw exception("Cannot move the iterator further");
			}

			setIterator++;

			if (setIterator == currentSetEnd)
			{
				mapIterator++;

				if (mapIterator != mapEnd)
				{
					setIterator = mapIterator->second.begin();
					currentSetEnd = mapIterator->second.end();
				}
			}

			return (*this);
		}

		Iterator operator++(int)
		{
			Iterator newIter(*this);
			operator++();
			return newIter;
		}

		NetworkAgentData* const& operator*()
		{
			if (setIterator == currentSetEnd)
			{
				return NULL;
			}

			return (*setIterator);
		}

		NetworkAgentData** operator->()
		{
			if (setIterator == currentSetEnd)
			{
				return NULL;
			}

			return (*setIterator);
		}

		Iterator() {}

		Iterator(const Iterator& src)
		{
			mapIterator = src.mapIterator;
			setIterator = src.setIterator;

			currentSetEnd = src.currentSetEnd;
			mapEnd = src.mapEnd;
		}

		Iterator& operator=(const Iterator& src)
		{
			if (this != &src)
			{
				this->mapIterator = src.mapIterator;
				this->setIterator = src.setIterator;

				this->currentSetEnd = src.currentSetEnd;
				this->mapEnd = src.mapEnd;
			}

			return (*this);
		}
	};

	Iterator BeginIterator()
	{
		Iterator begin;

		begin.mapIterator = networkAgents.begin();
		begin.mapEnd = networkAgents.end();

		if (begin.mapIterator != begin.mapEnd)
		{
			begin.setIterator = begin.mapIterator->second.begin();
			begin.currentSetEnd = begin.mapIterator->second.end();
		}

		return begin;
	}

	Iterator EndIterator()
	{
		Iterator end;

		end.mapIterator = networkAgents.end();
		end.mapEnd = networkAgents.end();

		return end;
	}

	NetworkAgentData* FindFirst(bool(*filter)(TAgent))
	{
		for (
			Iterator agentsIterator = BeginIterator();
			agentsIterator != EndIterator();
			agentsIterator++
			)
		{
			NetworkAgentData* candidate = const_cast<NetworkAgentData*>(*agentsIterator);
			if (filter(*(candidate->agent)))
			{
				return candidate;
			}
		}

		return NULL;
	}

#pragma region Functions for matrix operations

	/*
	 * The functions below are written to simplify
	 * operations with matrices. SwapColumns is much
	 * slower then SwapRows due to the how the matrix
	 * are stored (rowswisely).
	 *
	 * TODO:
	 * Addition and multiplication operations are
	 * not parallelled, though they should be.
	 */

	bool IsZero(double value)
	{
		return(abs(value) < 1e-7);
	}

	bool NonZero(double value)
	{
		return(!IsZero(value));
	}

	void SwapRows(vector<vector<double>>& matrix, int row1, int row2)
	{
		swap(matrix[row1], matrix[row2]);
	}

	void SwapColumns(vector<vector<double>>& matrix, int col1, int col2)
	{
		int fin = matrix.size();

		for (int i = 0; i < fin; i++)
		{
			swap(matrix[i][col1], matrix[i][col2]);
		}
	}

	void AddRow(vector<vector<double>>& matrix, int srcrow, int dstrow, double multiplier)
	{
		for (int i = 0; i < matrix[0].size(); i++)
		{
			matrix[dstrow][i] += matrix[srcrow][i] * multiplier;
		}
	}

	void AddColumn(vector<vector<double>>& matrix, int srccol, int dstcol, double multiplier)
	{
		for (int i = 0; i < matrix.size(); i++)
		{
			matrix[i][dstcol] += matrix[i][srccol] * multiplier;
		}
	}

	void MulRow(vector<vector<double>>& matrix, int row, double multiplier)
	{
		for (int i = 0; i < matrix[0].size(); i++)
		{
			matrix[row][i] *= multiplier;
		}
	}

	void MulColumn(vector<vector<double>>& matrix, int col, double multiplier)
	{
		for (int i = 0; i < matrix.size(); i++)
		{
			matrix[i][col] *= multiplier;
		}
	}

#pragma endregion

	/*
	 * bool SolveSLE(vector<vector<double>> matrix, vector<double> rhs) -
	 *     solves a system of linear equations.
	 *
	 *     Returns bool - true if the system if consistent, false otherwise
	 *
	 *     vector<vector<double>> matrix - the matrix of the system
	 *
	 *     vector<double> rhs - the right hand side column, if it's
	 *         a zero vector, the system is homogeneous
	 *
	 *     vector<double> solution - stores the solution of the SLE
	 *
	 *     The 'solution' values are unknown if the function
	 *     returns 'false'. If the system has an infinite amount
	 *     of solutions, the free variables will be set to 1.0
	 *     and the dependent ones will be calculated correspondingly
	 */
	bool SolveSLE(vector<vector<double>>& matrix, vector<double>& rhs, vector<double>& solution)
	{
		int equations_number = matrix.size();
		int vars_number = matrix[0].size();

		int current_row = 0;

		/*
		 * the loop below performs a downward
		 * walk of the gaussian elimination
		 */
		 // 'i' iterates through each variable
		for (int i = 0; i < vars_number; i++)
		{
			if (current_row == equations_number)
			{
				break;
			}

			/*
			 * first, we need to find a row, such that the
			 * i-th value in it (the coefficient of the current
			 * variable) is a non-zero value, we call this
			 * a nonzero_row
			 */
			int nonzero_row = current_row;
			while (IsZero(matrix[nonzero_row][i]))
			{
				nonzero_row++;
				// if there's no nonzero row, move to the next variable
				if (nonzero_row == equations_number)
				{
					break; // check out the 'if' block right below this while loop
				}
			}

			if (nonzero_row == equations_number)
			{
				continue; // jumps to the next iteration with i = i + 1
			}

			// replace the current row with the first found nonzero_row
			SwapRows(matrix, current_row, nonzero_row);
			swap(rhs[current_row], rhs[nonzero_row]);

			/*
			 * 'j' iterates through the rows below the current, nullifying
			 * corresponding values in the matrix
			 */
			for (int j = current_row + 1; j < equations_number; j++)
			{
				/*
				 * if the corresponding value is already nullified,
				 * move to the next row
				 */
				if (IsZero(matrix[j][i]))
				{
					continue;
				}

				// the number we multiply the current_row before subtracting it from the row j
				double coefficient = matrix[j][i] / matrix[current_row][i];

				for (int k = i; k < vars_number; k++)
				{
					matrix[j][k] -= matrix[current_row][k] * coefficient;
				}

				rhs[j] -= rhs[current_row] * coefficient;
			}

			current_row++;
		}

		if (current_row < equations_number)
		{
			/*
			 * in this case we need to check if
			 * any of the right-hand-side values
			 * below the current_row (including it)
			 * are non-zero - if this is true, then
			 * the SLE is inconsistent
			 */
			for (int i = current_row; i < equations_number; i++)
			{
				if (NonZero(rhs[i]))
				{
					return false;
				}
			}
		}

		// the only non-zero rows by now are 0 through current_now-1

		/*
		 * this int below holds
		 * vars_number-the number of calculated variables by the current moment
		 * made for convenience, so that we don't need to iterate
		 * until vars_number-the_other_possible_value (check out the inner loops)
		 */
		int calculated_offset = vars_number;

		/*
		 * this loop down below performs an upward walk
		 * of the gaussian elimination
		 * it handles both cases - when there's a unique
		 * solution and when there is an infinite amount
		 * of solutions - as mentioned in the description
		 * of the function, in this last case, it fills
		 * free variables with 1.0 values
		 */

		current_row--;

		for (; current_row >= 0; current_row--)
		{
			/*
			 * the leftmost nonzero value in any row
			 * for this moment can only be the diagonal
			 * value, which in-row index corresponds
			 * to the row number
			 */
			int leftmost_nonzero = current_row;
			while (IsZero(matrix[current_row][leftmost_nonzero]))
			{
				leftmost_nonzero++;
			}

			double lhs_calculated_part = 0.0;

			// in this loop we fill up free variables with 1.0 values
			for (int j = leftmost_nonzero + 1; j < calculated_offset; j++)
			{
				solution[j] = 1.0;
				lhs_calculated_part += matrix[current_row][j];
			}

			for (int j = calculated_offset; j < vars_number; j++)
			{
				lhs_calculated_part += matrix[current_row][j] * solution[j];
			}

			solution[leftmost_nonzero] = (rhs[current_row] - lhs_calculated_part) / matrix[current_row][leftmost_nonzero];

			/* okay, we've figured out the value for
			 * the current dependent variable
			 * (one for each row of the matrix)
			 * now we need to clean up the rows above
			 * the current one so that none of them
			 * contain this dependent variable as
			 * a summand (essentialy, this means that)
			 * we set the values in the corresponding
			 * column to 0, except the current row
			 */

			for (int j = current_row - 1; j >= 0; j--)
			{
				double multiplier = matrix[j][leftmost_nonzero] / matrix[current_row][leftmost_nonzero];
				for (int k = leftmost_nonzero; k < vars_number; k++)
				{
					matrix[j][k] -= matrix[current_row][k] * multiplier;
				}
			}

			calculated_offset = leftmost_nonzero;
		}

		/*
		 * this happens if we end up with a matrix
		 * with zero columns on the left by now
		 * (this is only possible if we were given
		 * them filled with zeros initially)
		 */
		for (int i = 0; i < calculated_offset; i++)
		{
			solution[i] = 1.0;
		}

		return true;
	}

public:
	void Add(TPosition position, TAgent* agent)
	{
		NetworkAgentData* newAgent = new NetworkAgentData(agent);

		for (
			Iterator agentsIterator = BeginIterator();
			agentsIterator != EndIterator();
			agentsIterator++
			)
		{
			NetworkAgentData* currentAgent = (*agentsIterator);

			if (currentAgent == newAgent)
			{
				continue;
			}

			DoNeighbor neighborChecker;

			if (neighborChecker(*(newAgent->agent), *(currentAgent->agent)))
			{
				newAgent->neighbors.insert(currentAgent);
				currentAgent->reverseNeighbors.insert(newAgent);
			}

			if (neighborChecker(*(currentAgent->agent), *(newAgent->agent)))
			{
				currentAgent->neighbors.insert(newAgent);
				newAgent->reverseNeighbors.insert(currentAgent);
			}
		}

		networkAgents[position].insert(newAgent);
		agentsCount++;
	}

	void Remove(TPosition position, TAgent* agent)
	{
		// if we do not have any agents at this position
		if (networkAgents.count(position) == 0)
		{
			return;
		}

		NetworkAgentData dummyToFind(agent);

		typename set<NetworkAgentData*, NetworkAgentDataComparator>::iterator agentData;
		agentData = networkAgents[position].find(&dummyToFind);

		// if there's no such agent
		if (agentData == networkAgents[position].end())
		{
			return;
		}

		NetworkAgentData* foundAgent = (*agentData);

		for (
			typename set<NetworkAgentData*>::iterator it = foundAgent->neighbors.begin();
			it != foundAgent->neighbors.end();
			it++
			)
		{
			(*it)->reverseNeighbors.erase(foundAgent);
		}

		for (
			typename set<NetworkAgentData*>::iterator it = foundAgent->reverseNeighbors.begin();
			it != foundAgent->reverseNeighbors.end();
			it++
			)
		{
			(*it)->neighbors.erase(foundAgent);
		}

		networkAgents[position].erase(foundAgent);
		delete foundAgent;
		agentsCount--;

		if (networkAgents[position].size() == 0)
		{
			networkAgents.erase(position);
		}
	}

	/*
	 * Calculates the diameter of the subgraph constructed of
	 * the vertices that satisfy the filter.
	 *
	 * When calling this function you need to be sure that
	 * the subgraph is connected, otherwise you're gonna
	 * get a diameter of any of its connectivity components.
	 */
	Measure GetDiameter(bool(*filter)(TAgent), bool(*findFirstFilter)(TAgent), Measure maxval)
	{
		Measure diameter = 0;

		vector<vector<int>> dp;

		dp.resize(agentsCount);
		for (int i = 0; i < agentsCount; i++)
		{
			dp[i].resize(agentsCount, maxval);
		}

		NetworkAgentData* root = FindFirst(findFirstFilter);

		if (root == NULL)
		{
			return ((Measure)0);
		}

		int subgraphsize = 0;
		map<NetworkAgentData*, int> agentsToIndicesMapping;

		queue<NetworkAgentData*> bfsQueue;

		agentsToIndicesMapping[root] = subgraphsize;
		dp[subgraphsize][subgraphsize] = 0;
		subgraphsize++;
		bfsQueue.push(root);

		while (!bfsQueue.empty())
		{
			NetworkAgentData* current = bfsQueue.front();
			bfsQueue.pop();

			for (
				typename set<NetworkAgentData*>::iterator neighbors = current->neighbors.begin();
				neighbors != current->neighbors.end();
				neighbors++
				)
			{
				NetworkAgentData* currentNeighbor = *neighbors;
				

				if (!filter(*(currentNeighbor->agent)))
				{
					continue;
				}

				// this condition is true if we've never met the agent before
				if (agentsToIndicesMapping.find(currentNeighbor) == agentsToIndicesMapping.end())
				{
					agentsToIndicesMapping[currentNeighbor] = subgraphsize;
					dp[subgraphsize][subgraphsize] = 0;
					subgraphsize++;
					bfsQueue.push(currentNeighbor);
				}
				// now we're sure that currentNeighbor has an index mapped to it

				CalculateDistance distanceCalculator;

				dp[agentsToIndicesMapping[current]][agentsToIndicesMapping[currentNeighbor]] =
					distanceCalculator(*(current->agent), *(currentNeighbor->agent));
			}
		}

		for (int currentAgent = 0; currentAgent < subgraphsize; currentAgent++)
		{
			for (int i = 0; i < subgraphsize; i++)
			{
				for (int j = 0; j < subgraphsize; j++)
				{
					dp[i][j] = min(dp[i][j], dp[i][currentAgent] + dp[currentAgent][j]);
				}
			}
		}

		for (int i = 0; i < subgraphsize; i++)
		{
			for (int j = 0; j < subgraphsize; j++)
			{
				diameter = max(diameter, dp[i][j]);
			}
		}

		return diameter;
	}

	vector<vector<int>> GetGeneratingFunctionValue(bool(*filter)(TAgent), bool(*findFirstFilter)(TAgent), Measure maxval)
	{
		// dp[i][j][k], i - ������ ���� � i, j - ����� � j, k - ����� ����� ����� k (���� ��������� ����� ���� = agentCount-1, ���� � ������� �� �������)
		vector<vector<vector<int>>> dp;
		vector<vector<int>> countVector;

		dp.resize(agentsCount);
		for (int i = 0; i < agentsCount; i++)
		{
			dp[i].resize(agentsCount);
			for (int j = 0; j < agentsCount; j++)
			{
				dp[i][j].resize(maxval, 0);
			}
		}

		NetworkAgentData* root = FindFirst(findFirstFilter);

		if (root == NULL)
		{
			// ��� ������ ������� ���� ����� ����� ������ �����
			// countVector[i][k], i - ����� ������, � ������ ���� ����������/�������������, k - ����� ����� ����� k
			// � countVector ��������� ��� ����, � �� ������ ���������� �����
			countVector.resize(1);	
			countVector[0].resize(1, 0);	
			return countVector;
		}

		int subgraphsize = 0;
		map<NetworkAgentData*, int> agentsToIndicesMapping;

		queue<NetworkAgentData*> bfsQueue;

		agentsToIndicesMapping[root] = subgraphsize;
		dp[subgraphsize][subgraphsize][0]++;
		subgraphsize++;
		bfsQueue.push(root);

		// Mapping agents with indices
		while (!bfsQueue.empty())
		{
			NetworkAgentData* current = bfsQueue.front();
			bfsQueue.pop();

			for (
				typename set<NetworkAgentData*>::iterator neighbors = current->neighbors.begin();
				neighbors != current->neighbors.end();
				neighbors++
				)
			{
				NetworkAgentData* currentNeighbor = *neighbors;

				if (!filter(*(currentNeighbor->agent)))
				{
					continue;
				}

				// this condition is true if we've never met the agent before
				if (agentsToIndicesMapping.find(currentNeighbor) == agentsToIndicesMapping.end())
				{
					agentsToIndicesMapping[currentNeighbor] = subgraphsize;
					dp[subgraphsize][subgraphsize][0]++;
					subgraphsize++;
					bfsQueue.push(currentNeighbor);
				}
				// now we're sure that currentNeighbor has an index mapped to it

				CalculateDistance distanceCalculator;
				dp[agentsToIndicesMapping[current]][agentsToIndicesMapping[currentNeighbor]][distanceCalculator(*(current->agent), *(currentNeighbor->agent))]++;
			}
		}

		// Calculate db
		for (int k = 1; k < maxval; k++)
		{
			for (
				typename std::map<NetworkAgentData*, int>::iterator startAgent = agentsToIndicesMapping.begin();
				startAgent != agentsToIndicesMapping.end();
				startAgent++
				)
			{
				for (
					typename std::map<NetworkAgentData*, int>::iterator finishAgent = agentsToIndicesMapping.begin();
					finishAgent != agentsToIndicesMapping.end();
					finishAgent++
					)
					{
						for (
						typename set<NetworkAgentData*>::iterator neighbor = finishAgent->first->neighbors.begin();
						neighbor != finishAgent->first->neighbors.end();
						neighbor++
						)
						{
							dp[startAgent->second][finishAgent->second][k] += dp[startAgent->second][agentsToIndicesMapping[*neighbor]][k - 1];
						}
					}
			}
		}

		// Init countVector
		// ��� ������ ������� ���� ����� ����� ������ �����
		// countVector[i][k], i - ����� ������, � ������ ���� ����������/�������������, k - ����� ����� ����� k
		// � countVector ��������� ��� ����, � �� ������ ���������� �����
			
		countVector.resize(agentsCount);
		for (int i = 0; i < agentsCount; i++)
		{
			countVector[i].resize(maxval, 0);
		}

		// Calculate countVector
		for (int i = 0; i < agentsCount; i++)
		{
			for (int k = 0; k < maxval; k++)
			{
				for (int j = 0; j < agentsCount; j++)
				{
					countVector[i][k] += dp[i][j][k];
				}

			}	
		}

		return countVector;
	}

	vector<pair<TAgent*, double>> GetPageRank(bool(*filter)(TAgent), bool(*findFirstFilter)(TAgent))
	{
		vector<pair<TAgent*, double>> result;

		vector<vector<double>> ranks;

		ranks.resize(agentsCount);

		for (int i = 0; i < agentsCount; i++)
		{
			ranks[i].resize(agentsCount, 0.0);
		}

		NetworkAgentData* root = FindFirst(findFirstFilter);

		if (root == NULL)
		{
			throw "No agents correspond to the filter";
		}

		int subgraphsize = 0;
		map<NetworkAgentData*, int> agentsToIndicesMapping;

		queue<NetworkAgentData*> bfsQueue;

		agentsToIndicesMapping[root] = subgraphsize;
		subgraphsize++;
		bfsQueue.push(root);

		vector<int> filteredNeighbors;
		result.push_back(make_pair(root->agent, 0.0));

		while (!bfsQueue.empty())
		{
			NetworkAgentData* current = bfsQueue.front();
			bfsQueue.pop();

			int currentFilteredNeighbors = 0;

			for (
				typename set<NetworkAgentData*>::iterator neighbors = current->neighbors.begin();
				neighbors != current->neighbors.end();
				neighbors++
				)
			{
				NetworkAgentData* currentNeighbor = *neighbors;

				if (!filter(*(currentNeighbor->agent)))
				{
					continue;
				}

				currentFilteredNeighbors++;

				// this condition is true if we've never met the agent before
				if (agentsToIndicesMapping.find(currentNeighbor) == agentsToIndicesMapping.end())
				{
					agentsToIndicesMapping[currentNeighbor] = subgraphsize;
					subgraphsize++;
					bfsQueue.push(currentNeighbor);
					result.push_back(make_pair(currentNeighbor->agent, 0.0));
				}
				// now we're sure that currentNeighbor has an index mapped to it

				// since these guys are neighbors, we set the value to 1.0
				ranks[agentsToIndicesMapping[currentNeighbor]][agentsToIndicesMapping[current]] = 1.0;
			}

			filteredNeighbors.push_back(currentFilteredNeighbors);
		}

		if (subgraphsize == 1)
		{
			result.push_back(make_pair(root->agent, 0.0));
			return result;
		}

		//vector<int> inversions;

		for (int i = 0; i < subgraphsize; i++)
		{
			//inversions.push_back(i);

			ranks[i][i] = -1.0;

			for (int j = 0; j < subgraphsize; j++)
			{
				if (i != j)
				{
					ranks[j][i] /= ((double)filteredNeighbors[i]);
				}
			}
		}

		//reverse(ranks.begin(), ranks.begin()+subgraphsize);

		result.resize(subgraphsize);
		for (
			typename std::map<NetworkAgentData*, int>::iterator it = agentsToIndicesMapping.begin();
			it != agentsToIndicesMapping.end();
			it++
			)
		{
			result[it->second].first = it->first->agent;
		}

		vector<double> rightHandSide;
		rightHandSide.resize(subgraphsize, 0.0);

		vector<double> sleSolution;
		sleSolution.resize(subgraphsize);

		// now we have a matrix built

		if (!SolveSLE(ranks, rightHandSide, sleSolution))
		{
			throw exception("The SLE for calculating the PageRank is invalid");
		}

		for (int i = 0; i < subgraphsize; i++)
		{
			result[i].second = sleSolution[i];
		}

		return result;
	}
	/*
	void Test()
	{
		Iterator it = BeginIterator();

		while( it != EndIterator() )
			cout << it->neighbors.size() << endl;
	}
	*/
	
	/*
	vector<pair<TAgent*, double>> FBetwCentr(bool(*filter)(TAgent), bool(*findFirstFilter)(TAgent))
	{
		vector<pair<TAgent*, double>> cb;
		map<int, NetworkAgentData*> agentsToIndicesMapping;
		NetworkAgentData* currentNeighbor;
		NetworkAgentData* temp;

		NetworkAgentData* root = FindFirst(findFirstFilter);
		if (root == NULL)
		{
			throw "No agents correspond to the filter";
		}

		int numbVrtx = 0;
		int v;
		int w;

		//agentsToIndicesMapping[root] = numbVrtx;
		agentsToIndicesMapping[numbVrtx] = root;


		cb.push_back(make_pair(root->agent, 0.0));


		stack<int> S;
		vector<int> d;
		vector<int> prev;
		queue<int> Q;
		vector<float> sigma;
		map< int, vector<int> > P; // map for mapping neighbours of a vertex to itself



		// Initialize various BFS variables
		for (int i = 0; i < agentsCount; i++)
		{
			d.push_back(-1);
			sigma.push_back(0);
			prev.push_back(0);
		}

		sigma[numbVrtx] = 1;
		d[numbVrtx] = 0;
		Q.push(numbVrtx);
		

		// Perform the BFS
		while (!Q.empty())
		{
			
			v = Q.front();
			Q.pop();
			S.push(v);
			temp = agentsToIndicesMapping[v] ;

		

			// Traverse the neighbours of v
			for (typename set<NetworkAgentData*>::iterator neighbors = temp->neighbors.begin();
				neighbors != temp->neighbors.end();
				neighbors++)
			{
				NetworkAgentData* currentNeighbor = *neighbors;
				
				
				if (agentsToIndicesMapping.find(numbVrtx + 1) == agentsToIndicesMapping.end()) //тут должно быть условие если d[v]<0
				{
					agentsToIndicesMapping[numbVrtx + 1] = currentNeighbor;
					w = numbVrtx ++;

					if (d[w] < 0)
					{

						
						Q.push(w);
						d[w] = d[v] + 1;

					}
				}

				// Update sigma value of w if path through v was shortest
				if (d[w] == d[v] + 1)
				{
					sigma[w] = sigma[w] + sigma[v];
					P[w].push_back(v);
				}
			}
		}

		vector<float> delta;
		for (int i = 0; i < agentsCount; i++)
		{
			delta.push_back(0);
		}


		// Pop elements out of the stack, starting from terminal node 
		// work backward frontier b frontier, computing delta values
		while (!S.empty()) {
			int w = S.top();
			S.pop();

			vector<int>::iterator it; // we needed an iterator as we used an STL vector
			for (it = P[w].begin(); it != P[w].end(); it++) {
				int v = *it;
				delta[v] += ((sigma[v] / sigma[w]) * (1 + delta[w]));

			}
			if (w != 0) {
				cb[w].second += delta[w] / 2;
			}
		}



		for (
			typename std::map<int, NetworkAgentData*>::iterator it = agentsToIndicesMapping.begin();
			it != agentsToIndicesMapping.end();
			it++
			)

		{
			cb[it->first].first = it->second->agent;
		}

		return cb;
	}
	
	vector<pair<TAgent*, double>> FindBetwCentr(bool(*filter)(TAgent), bool(*findFirstFilter)(TAgent))
	{
		vector<pair<TAgent*, double>> result;
		NetworkAgentData* root = FindFirst(findFirstFilter);

		int numbVrtx = 0; //номер рассматриваемой вершины

		map<NetworkAgentData*, int> agentsToIndicesMapping;
		queue<NetworkAgentData*> Q;
		stack <NetworkAgentData*> S; //может тут тоже NetworkAgentData* ?


		vector <long> d;						// A vector storing shortest distance estimates
		vector <long> sigma;				// sigma is the number of shortest paths
		vector <double> delta;							// A vector storing dependency of the source vertex on all other vertices
		vector< vector <unsigned long> > Pred;		// A list of predecessors of all vertices 

		NetworkAgentData* u;
		NetworkAgentData* currentNeighbor;

		vector<int> filteredNeighbors;
		int currentFilteredNeighbors;

		d.assign(agentsCount, ULONG_MAX);
		Pred.assign(agentsCount, vector <unsigned long>(0, 0));
		sigma.assign(agentsCount, 0);
		delta.assign(agentsCount, 0);

		if (root == NULL)
		{
			throw "No agents correspond to the filter";
		}

		agentsToIndicesMapping[root] = numbVrtx; // root - это ключ, по которому доступ к mapping

		d[numbVrtx] = 0;
		sigma[numbVrtx] = 1;
													// нашли первого, удовл-его фильтру, сказали, что subgraph - будет его номер
		//numbVrtx++;   //?

		Q.push(root);	// в очередь добавили 1го агента, удовл-го фильтру.

		result.push_back(make_pair(root->agent, 0.0));	
		

		
		

		while (!Q.empty())
		{
			u = Q.front(); //1ый элемент из очереди = u
			d[numbVrtx] = 0;
			sigma[agentsToIndicesMapping[u]] = 1;
			Q.pop();
			S.push(u);

			currentFilteredNeighbors = 0;

			for (
				typename set<NetworkAgentData*>::iterator neighbors = u->neighbors.begin();
				neighbors != u->neighbors.end();
				neighbors++
				)
			{
				currentNeighbor = *neighbors;


				currentFilteredNeighbors++;


				// this condition is true if we've never met the agent before
				if (agentsToIndicesMapping.find(currentNeighbor) == agentsToIndicesMapping.end()) //тут должно быть условие если d[v]<0
				{
					numbVrtx++;
					agentsToIndicesMapping[currentNeighbor] = numbVrtx; //subgraphsize - это же номер соседа?
					
					Q.push(currentNeighbor); //добавляем соседа в Q

					d[numbVrtx] = d[agentsToIndicesMapping[u]] + 1; // расстояние

					result.push_back(make_pair(currentNeighbor->agent, 0.0)); //это хз
				}

				if (d[numbVrtx] = d[agentsToIndicesMapping[u]] + 1)
				{
					sigma[agentsToIndicesMapping[currentNeighbor]] += sigma[agentsToIndicesMapping[u]];
					Pred[agentsToIndicesMapping[currentNeighbor]].push_back(agentsToIndicesMapping[u]);
				}


			}
		}

			
		
			//result.resize(numbVrtx);
			for (
				typename std::map<NetworkAgentData*, int>::iterator it = agentsToIndicesMapping.begin();
				it != agentsToIndicesMapping.end();
				it++
				)

			{
				result[it->second].first = it->first->agent;
			}



			while(!S.empty())
			{
				u = S.top();
				S.pop();

				for (int j = 0; j < Pred[agentsToIndicesMapping[u]].size(); j++)
				{
					delta[Pred[agentsToIndicesMapping[u]][j]] += (sigma[Pred[agentsToIndicesMapping[u]][j]] / sigma[agentsToIndicesMapping[u]]) * (1 + delta[agentsToIndicesMapping[u]]);
				}


				//if (agentsToIndicesMapping[root] != agentsToIndicesMapping[u])
				
				result[agentsToIndicesMapping[u]].second += delta[agentsToIndicesMapping[u]];
				
			}

			return result;
	}
			
*/

	vector<pair<TAgent*, double>> BetwCentral(bool(*filter)(TAgent))
	{
		vector<pair<TAgent*, double>> result;
		NetworkAgentData* root = FindFirst(filter);// root - это ключ, по которому доступ к mapping

		int numbVrtx = 0; //номер рассматриваемой вершины
		

		map< NetworkAgentData*, int> agentsToIndicesMapping;
		map< int, NetworkAgentData*> indicesToAgentsMapping;
		queue<NetworkAgentData*> Q;
		stack <NetworkAgentData*> S; //может тут тоже NetworkAgentData* ?
		vector <double> CB;

		vector <long> d;						// A vector storing shortest distance estimates
		vector <long> sigma;				// sigma is the number of shortest paths
		vector <double> delta;							// A vector storing dependency of the source vertex on all other vertices
		vector< vector <unsigned long> > Pred;		// A list of predecessors of all vertices 

		NetworkAgentData* u;
		NetworkAgentData* currentNeighbor;

		vector<int> filteredNeighbors;
	
		agentsToIndicesMapping[root] = numbVrtx;
		indicesToAgentsMapping[numbVrtx] = root;

		if (root == NULL)
		{
			throw "No agents correspond to the filter";
		}
		Q.push(root);

		while (!Q.empty())
		{
			u = Q.front(); //1ый элемент из очереди = u
			Q.pop();

			for (
				typename set<NetworkAgentData*>::iterator neighbors = u->neighbors.begin();
				neighbors != u->neighbors.end();
				neighbors++
				)
			{
				currentNeighbor = *neighbors;
				
				// this condition is true if we've never met the agent before
				if (agentsToIndicesMapping.find(currentNeighbor) == agentsToIndicesMapping.end()) //тут должно быть условие если d[v]<0
				{
					numbVrtx++;
					agentsToIndicesMapping[currentNeighbor] = numbVrtx;
					indicesToAgentsMapping[numbVrtx] = currentNeighbor;

					Q.push(currentNeighbor); //добавляем соседа в Q
				}
			}
		}
		
		CB.assign(agentsToIndicesMapping.size(), 0);

		queue<int> Q2;

		for (int i = 0; i < agentsToIndicesMapping.size(); i++)
		{
			d.assign(agentsToIndicesMapping.size(), ULONG_MAX);
			Pred.assign(agentsToIndicesMapping.size(), vector <unsigned long>(0, 0));
			sigma.assign(agentsToIndicesMapping.size(), 0);
			delta.assign(agentsToIndicesMapping.size(), 0);
			root = indicesToAgentsMapping[i];
			
			d[i] = 0;
			sigma[i] = 1;
		
			Q2.push(i);

			

			while (!Q2.empty())
			{
				u = indicesToAgentsMapping[Q2.front()]; //1ый элемент из очереди = u
				Q2.pop();
				S.push(u);

				for (
					typename set<NetworkAgentData*>::iterator neighbors = u->neighbors.begin();
					neighbors != u->neighbors.end();
					neighbors++
					)
				{
					currentNeighbor = *neighbors;

					// this condition is true if we've never met the agent before
					if (d[agentsToIndicesMapping[currentNeighbor]] == ULONG_MAX) //тут должно быть условие если d[v]<0
					{
						d[agentsToIndicesMapping[currentNeighbor]] = d[agentsToIndicesMapping[u]] + 1; // расстояние
						Q2.push(agentsToIndicesMapping[currentNeighbor]); //добавляем соседа в Q
					}

					if (d[agentsToIndicesMapping[currentNeighbor]] == d[agentsToIndicesMapping[u]] + 1)
					{
						sigma[agentsToIndicesMapping[currentNeighbor]] += sigma[agentsToIndicesMapping[u]];
						Pred[agentsToIndicesMapping[currentNeighbor]].push_back(agentsToIndicesMapping[u]);
					}
				}
			}


			while (!S.empty())
			{
				u = S.top();
				S.pop();

				for (int j = 0; j < Pred[agentsToIndicesMapping[u]].size(); j++)
				{
					delta[Pred[agentsToIndicesMapping[u]][j]] += ((double)sigma[Pred[agentsToIndicesMapping[u]][j]] / sigma[agentsToIndicesMapping[u]]) * (1 + delta[agentsToIndicesMapping[u]]);
				}

				if (agentsToIndicesMapping[u] != i)
					CB[agentsToIndicesMapping[u]] += delta[agentsToIndicesMapping[u]];
			}

			// Clear data for the next run
			Pred.clear();
			d.clear();
			sigma.clear();
			delta.clear();

		} // End For 
		

		for (
			typename std::map<NetworkAgentData*, int>::iterator it = agentsToIndicesMapping.begin();
			it != agentsToIndicesMapping.end();
			it++
			)

		{
			result[it->second].first = it->first->agent;
		}

		//result.push_back(make_pair(root->agent, 0.0));
		return result;
	}

	vector< vector <int> > GreedyAl(bool(*filter)(TAgent))
	{
		vector< vector <int>> result;
		NetworkAgentData* root = FindFirst(filter);// root - это ключ, по которому доступ к mapping

		int numbVrtx = 0; //номер рассматриваемой вершины

		int max[3] = { 0, 0, 0 };
		map< NetworkAgentData*, int> agentsToIndicesMapping;
		map< int, NetworkAgentData*> indicesToAgentsMapping;
		queue<NetworkAgentData*> Q;
		vector<vector<int>> dp;
//		CalculateDistance distanceCalculator;
		
		NetworkAgentData* u;
		NetworkAgentData* currentNeighbor;

		agentsToIndicesMapping[root] = numbVrtx;
		indicesToAgentsMapping[numbVrtx] = root;

		if (root == NULL)
		{
			throw "No agents correspond to the filter";
		}
		Q.push(root);

		while (!Q.empty())
		{
			u = Q.front(); //1ый элемент из очереди = u
			Q.pop();

			for (
				typename set<NetworkAgentData*>::iterator neighbors = u->neighbors.begin();
				neighbors != u->neighbors.end();
				neighbors++
				)
			{
				currentNeighbor = *neighbors;

				// this condition is true if we've never met the agent before
				if (agentsToIndicesMapping.find(currentNeighbor) == agentsToIndicesMapping.end()) //тут должно быть условие если d[v]<0
				{
					numbVrtx++;
					agentsToIndicesMapping[currentNeighbor] = numbVrtx;
					indicesToAgentsMapping[numbVrtx] = currentNeighbor;

					Q.push(currentNeighbor); //добавляем соседа в Q
				}
				CalculateDistance distanceCalculator;

				dp[agentsToIndicesMapping[current]][agentsToIndicesMapping[currentNeighbor]] =
					distanceCalculator(*(current->agent), *(currentNeighbor->agent));
			}
		}

		dp.resize(agentsToIndicesMapping.size());
		for (int i = 0; i < agentsToIndicesMapping.size(); i++)
		{
			dp[i].resize(agentsToIndicesMapping.size(), 0);
		}

		for (int i = 0; i < agentsToIndicesMapping.size(); i++)
		{
			for (int j = 0; j < agentsToIndicesMapping.size(); j++)
			{
				dp[i][j] = distanceCalculator(*(indicesToAgentsMapping[i]->agent), *(indicesToAgentsMapping[j]->agent));
				if (dp[i][j] > max[0])
				{
					max[0] = dp[i][j];
					max[1] = i;
					max[2] = j;
					
				}
			}
		}

		for (int currentAgent = 0; currentAgent < agentsToIndicesMapping.size(); currentAgent++)
		{
			for (int i = 0; i < agentsToIndicesMapping.size(); i++)
			{
				for (int j = 0; j < agentsToIndicesMapping.size(); j++)
				{
					dp[i][j] = min(dp[i][j], dp[i][currentAgent] + dp[currentAgent][j]);
				}
			}
		}

		return result;
	}

	Network() {}
};

#endif /* MANET_PROJECT_NETWORK_HEADER_DEFINED */