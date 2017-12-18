#include "SimModel.h"
#include<fstream>
#include<cassert>
#include<iomanip>
#include "element.h"


SimModel::SimModel(std::string simulationName)
{
	simModelName = simulationName;			
}

bool SimModel::createSimModelFromInp(std::string strInpFileIn)
{
	std::ifstream myInpFile;
	std::string sInpFile;
	myInpFile.open(strInpFileIn);
	myInpFile.seekg(0, myInpFile.end);

	sInpFile.resize(myInpFile.tellg());
	myInpFile.seekg(0, std::ios::beg);
	myInpFile.read(&sInpFile[0], sInpFile.size());

	

	

	//--------------------------
	//Crude method-1
	//Create - Node Array
	//Create - Element Array
	//Read element type
	//--------------------------

	//------------
	//Find Nodes
	//------------
	std::smatch m;
	std::regex e("^\\*node.*?\\n([\\s\\S]*?)\\n\\*", std::regex_constants::icase);
	std::string s = sInpFile;
	bool isMatch = std::regex_search(s, m, e);
	std::string sNodeData;

	if (isMatch) {
		sNodeData = m[1];
	}
	else {
		std::cout << "No Node data found.";
		return false;
	}


	//------------
	//Find Elements
	//------------
	std::string sElementType;
	e = std::regex("^\\*element.*?,.*?type=(.*)?\\n([\\s\\S]*?)\\n\\*", std::regex_constants::icase);
	isMatch = std::regex_search(s, m, e);
	std::string sElementData;
	if (isMatch) {
		sElementType = m[1];
		sElementData = m[2];
		//std::cout << sElementType << std::endl<<sElementData;
	}
	else {
		std::cout << "No Element data found. The program will now exit";
		return false;
	}

	

	//---------------
	//Get the element object from the element type - not implemented
	//---------------
	int numNodesPerElement;
	modelData.numNodes = getNumNodes(sNodeData);
	modelData.numElements = getNumElements(sElementData, numNodesPerElement);
	modelData.iDimension = 2;

	if (modelData.elementLibrary.find(sElementType) == modelData.elementLibrary.end()) {
		std::cout << "Element type not supported";
		return false;
	}
	else {

		modelData.currElement = modelData.elementLibrary[sElementType];
	}
	


	//gsl_matrix * nodeCoordArray = gsl_matrix_alloc(numNodes, 2); //Allocate Node array
	modelData.nodeCoordArray = gsl_matrix_alloc(modelData.numNodes, 2); //Allocate Node array
	modelData.elementConnectivityArray = gsl_matrix_int_alloc(modelData.numElements, numNodesPerElement); //Allocate Node array

	bool isNodeArraySuccess = populateNodeArray(sNodeData, modelData.nodeCoordArray, modelData.iDimension);
	
	if (!isNodeArraySuccess)
	{
		std::cout << "Failed to create array of node coords";
		return false;
	}
	else std::cout << "Processed *node block";

	bool isElementArraySuccess = populateELementArray(sElementData, modelData.elementConnectivityArray, modelData.numElements, numNodesPerElement);
	if (!isElementArraySuccess)
	{
		std::cout << "Failed to create array of elements";
		return false;
	}
	else std::cout << "Processed *element block";


	//-------------------
	//Create Node and Element set map
	//-------------------

	//------------
	//Find Elements sets
	//------------	

	try {
		std::regex elsetRegex("(^\\*elset.*\\n[\\s\\S]*?)(?=\\n\\*)", std::regex_constants::icase);
		//std::regex elsetRegex("^\\*elset", std::regex_constants::icase);

		std::sregex_iterator next(s.begin(), s.end(), elsetRegex);
		std::sregex_iterator end;
		while (next != end) {
			m = *next;
			std::cout << m.str() << "\n";
			modelData.addElementSetToMap(m[1]);
			next++;
		}
	}
	catch (std::regex_error& err) {
		// Syntax error in the regular expression
		std::cout << "No element sets detected";
	}

	//------------
	//Find Node sets
	//------------
	try {
		std::regex elsetRegex("(^\\*nset.*\\n[\\s\\S]*?)(?=\\n\\*)", std::regex_constants::icase);
		std::sregex_iterator next(s.begin(), s.end(), elsetRegex);
		std::sregex_iterator end;
		while (next != end) {
			m = *next;
			std::cout << m.str() << "\n";
			modelData.addNodeSetToMap(m[1]);
			next++;
		}
	}
	catch (std::regex_error& err) {
		// Syntax error in the regular expression
		std::cout << "No node sets detected";
	}
	

	//---------------------
	//Find amplitudes
	//--------------------

	try {
		std::regex elsetRegex("(^\\*amplitude.*\\n[\\s\\S]*?)(?=\\n\\*)", std::regex_constants::icase);
		std::sregex_iterator next(s.begin(), s.end(), elsetRegex);
		std::sregex_iterator end;
		while (next != end) {
			m = *next;
			std::cout << m.str() << "\n";
			modelData.addAmplitudesToMap(m[1]);
			next++;
		}
	}
	catch (std::regex_error& err) {
		// Syntax error in the regular expression
		std::cout << "No node sets detected";
	}

	std::cout << "Model Created\n";

	//--------------------------------------------------------------
	//History Data
	//--------------------------------------------------------------
	

	//-----------------
	//Process step data
	//-----------------

	try {
		std::regex elsetRegex("(^\\*step.*\\n[\\s\\S]*?)(?=\\n\\*end step)", std::regex_constants::icase);
		std::sregex_iterator next(s.begin(), s.end(), elsetRegex);
		std::sregex_iterator end;
		while (next != end) {
			m = *next;
			std::cout << m.str() << "\n";
			historyData.addSteps(m[1]);
			next++;
		}
	}
	catch (std::regex_error& err) {
		// Syntax error in the regular expression
		std::cout << "No Step data found";
	}

	std::cout << "History data created";


}

int SimModel::getNumNodes(const std::string & nodeBlockString)
{
	int numNodes = 0;
	for (std::string::size_type i = 0; i < nodeBlockString.size(); i++)
	{
		if (nodeBlockString[i] == '\n') numNodes++;
	}
	numNodes++;
	if (numNodes > 1) return numNodes;
	return -1;
}

int SimModel::getNumElements(const std::string & elementBlockString, int & numNodesPerElement)
{
	int numELements = 0;
	numNodesPerElement = 0;
	for (std::string::size_type i = 0; i < elementBlockString.size(); i++)
	{
		if (elementBlockString[i] == ',') numNodesPerElement++;
		if (elementBlockString[i] == '\n') break;
	}

	for (std::string::size_type i = 0; i < elementBlockString.size(); i++)
	{
		if (elementBlockString[i] == '\n') numELements++;
	}
	if (numNodesPerElement > 0)  numELements++;

	if (numELements > 0) return numELements;
	return -1;
}

bool SimModel::populateNodeArray(const std::string & nodeBlockString, gsl_matrix * nodeArray, int iDimension)
{
	std::string nodeLine;
	std::string nodeTemp;
	double temp;
	int cnt = 0;
	int nodeCnt = 0;


	if (iDimension == 2)
	{
		std::istringstream nodeData(nodeBlockString);
		while (std::getline(nodeData, nodeLine, '\n'))
		{
			std::istringstream nodeLineStream(nodeLine);
			cnt = 0;
			while (std::getline(nodeLineStream, nodeTemp, ','))
			{
				temp = std::stod(nodeTemp);
				if (cnt == 1) gsl_matrix_set(nodeArray, nodeCnt, 0, temp);
				if (cnt == 2) gsl_matrix_set(nodeArray, nodeCnt, 1, temp);
				cnt++;
			}
			nodeCnt++;
		}

	}
	else if (iDimension == 3)
	{
		std::istringstream nodeData(nodeBlockString);
		while (std::getline(nodeData, nodeLine, '\n'))
		{
			std::istringstream nodeLineStream(nodeLine);
			cnt = 0;
			while (std::getline(nodeLineStream, nodeTemp, ','))
			{
				temp = std::stod(nodeTemp);
				if (cnt == 1) gsl_matrix_set(nodeArray, nodeCnt, 0, temp);
				if (cnt == 2) gsl_matrix_set(nodeArray, nodeCnt, 1, temp);
				if (cnt == 3) gsl_matrix_set(nodeArray, nodeCnt, 2, temp);
				cnt++;
			}
			nodeCnt++;
		}

	}


	//Add checks for validity
	return true;

	return false;
}

bool SimModel::populateELementArray(const std::string & elementBlockString, gsl_matrix_int * elementArray, const int numElements, const int numNodesPerElement)
{
	std::string elementLine;
	std::string elementTemp;
	int temp;
	int cnt = 0;
	int elemCnt = 0;

	std::istringstream elemData(elementBlockString);
	while (std::getline(elemData, elementLine, '\n'))
	{
		std::istringstream elemLineStream(elementLine);
		cnt = 0;
		while (std::getline(elemLineStream, elementTemp, ','))
		{
			temp = std::stoi(elementTemp);
			if (cnt>0)gsl_matrix_int_set(elementArray, elemCnt, cnt - 1, temp);
			cnt++;
		}
		elemCnt++;
	}

	//Add checks for validity
	return true;



	return false;
}

SimModel::~SimModel()
{
}










//--------------------------------------------------------------------
//Model Data
//--------------------------------------------------------------------

ModelData::ModelData()
{

	//Create the element library
	//Element anyElement;	
	elementLibrary["CPE8"] =new Element_CPE8();	
	elementLibrary["CPE6"] = new Element_CPE6();
	elementLibrary["CPE3"] = new Element_CPE3();
}

bool ModelData::addElementSetToMap(std::string elsetStr)
{
	//process line-1
	std::string elsetLine;
	std::string elTemp, elIdTemp;
	int numParams;
	//std:string elset;
	int cnt = 0;
	std::string sElsetName;
	bool isGenerateOn = false;
	int genStart, genEnd, genIter;
	bool isMatch;

	std::vector<int> tempElsetVec;

	std::smatch m;
	std::regex e("elset=(.*?),", std::regex_constants::icase);
	std::string s = "";
	gsl_vector_int * elSetArray;
	int elmCnt=0;


	std::istringstream elementSetData(elsetStr);
	while (std::getline(elementSetData, elsetLine, '\n'))
	{

		cnt++;
		if (cnt == 1)
		{
			//Get elset name
			e = std::regex("elset=(.*?),", std::regex_constants::icase);
			s = elsetLine;
			isMatch = std::regex_search(s, m, e);
			if (isMatch) sElsetName = m[1];			
			else
			{
				e = std::regex("elset=(.*)", std::regex_constants::icase);
				s = elsetLine;
				isMatch = std::regex_search(s, m, e);
				if (isMatch) sElsetName = m[1];

				else {
					std::cout << "Element set creation failed";
					return false;
				}
					
			}

			//Check whether generate is turned on
			e = std::regex("generate", std::regex_constants::icase);
			s = elsetLine;
			isMatch = std::regex_search(s, m, e);
			if (isMatch) isGenerateOn = true;

		}
		else
		{
			std::istringstream elsetLineStream(elsetLine);

			if (isGenerateOn)
			{
				if (cnt == 2)
				{
					int cnt1 = 0;
					while (std::getline(elsetLineStream, elIdTemp, ','))
					{
						cnt1++;
						int temp = std::stoi(elIdTemp);
						if (cnt1 == 1) genStart = temp;
						if (cnt1 == 2) genEnd = temp;
						if (cnt1 == 3) genIter = temp;
					}
					if (cnt1 < 3) {
						std::cout << "Element set creation failed";
						return false;
					}
					int cnt2 = 0;
					for (int i = genStart; i <= genEnd; i = i + genIter) {
						cnt2++;
					}
					elSetArray = gsl_vector_int_alloc(cnt2);
					cnt2 = 0;
					for (int i = genStart; i <= genEnd; i = i + genIter) {
						gsl_vector_int_set(elSetArray, cnt2, i);
						cnt2++;
					}
					elementSetMap[sElsetName] = elSetArray;
					return true;
				}
			}
			else
			{
				while (std::getline(elsetLineStream, elIdTemp, ','))
				{
					tempElsetVec.push_back(std::stoi(elIdTemp));
				}
			}
		}
	}

	if (!isGenerateOn)
	{
		elSetArray = gsl_vector_int_alloc(tempElsetVec.size());
		std::vector<int>::const_iterator iter;
		cnt = 0;
		for (iter = tempElsetVec.begin(); iter != tempElsetVec.end(); ++iter)
		{
			gsl_vector_int_set(elSetArray, cnt, *iter);
			cnt++;				
		}
		elementSetMap[sElsetName] = elSetArray;
		return true;
	}

	std::cout << "Element set creation failed";
	return false;
}

bool ModelData::addNodeSetToMap(std::string nodeSetStr)
{
	//process line-1
	std::string nodesetLine;
	std::string nodeTemp, nodeIdTemp;
	int numParams;
	//std:string nodeset;
	int cnt = 0;
	std::string sNodeSetName;
	bool isGenerateOn = false;
	int genStart, genEnd, genIter;
	bool isMatch;

	std::vector<int> tempnodesetVec;

	std::smatch m;
	std::regex e("nset=(.*?),", std::regex_constants::icase);
	std::string s = "";
	gsl_vector_int * nodeSetArray;
	int nodeCnt = 0;


	std::istringstream nodeSetData(nodeSetStr);
	while (std::getline(nodeSetData, nodesetLine, '\n'))
	{

		cnt++;
		if (cnt == 1)
		{
			//Get nodeset name
			e = std::regex("nset=(.*?),", std::regex_constants::icase);
			s = nodesetLine;
			isMatch = std::regex_search(s, m, e);
			if (isMatch) sNodeSetName = m[1];
			else
			{
				e = std::regex("nset=(.*)", std::regex_constants::icase);
				s = nodesetLine;
				isMatch = std::regex_search(s, m, e);
				if (isMatch) sNodeSetName = m[1];
				else
				{
					std::cout << "Node set creation failed";
					return false;
				}
			}

			//Check whether generate is turned on
			e = std::regex("generate", std::regex_constants::icase);
			s = nodesetLine;
			isMatch = std::regex_search(s, m, e);
			if (isMatch) isGenerateOn = true;

		}
		else
		{
			std::istringstream nodesetLineStream(nodesetLine);

			if (isGenerateOn)
			{
				if (cnt == 2)
				{
					int cnt1 = 0;
					while (std::getline(nodesetLineStream, nodeIdTemp, ','))
					{
						cnt1++;
						int temp = std::stoi(nodeIdTemp);
						if (cnt1 == 1) genStart = temp;
						if (cnt1 == 2) genEnd = temp;
						if (cnt1 == 3) genIter = temp;
					}		   
					if (cnt1 < 3) {
						std::cout << "Node set creation failed";
						return false;
					}

					int cnt2 = 0;
					for (int i = genStart; i <= genEnd; i = i + genIter) {
						cnt2++;
					}
					nodeSetArray = gsl_vector_int_alloc(cnt2);
					cnt2 = 0;
					for (int i = genStart; i <= genEnd; i = i + genIter) {
						gsl_vector_int_set(nodeSetArray, cnt2, i);
						cnt2++;
					}
					nodeSetMap[sNodeSetName] = nodeSetArray;
					return true;
				}
			}
			else
			{
				while (std::getline(nodesetLineStream, nodeIdTemp, ','))
				{
					tempnodesetVec.push_back(std::stoi(nodeIdTemp));
				}
			}
		}
	}

	if (!isGenerateOn)
	{
		nodeSetArray = gsl_vector_int_alloc(tempnodesetVec.size());
		std::vector<int>::const_iterator iter;
		cnt = 0;
		for (iter = tempnodesetVec.begin(); iter != tempnodesetVec.end(); ++iter)
		{
			gsl_vector_int_set(nodeSetArray, cnt, *iter);
			cnt++;
		}
		nodeSetMap[sNodeSetName] = nodeSetArray;
		return true;
	}

	std::cout << "Node set creation failed";
	return false;

}

bool ModelData::addAmplitudesToMap(std::string ampStringIn)
{
	//process line-1
	std::string ampLine;
	std::string ampTimeTemp, ampDataTemp;
	int numParams;
	//std:string nodeset;
	int cnt = 0;
	std::string sAmpName;
	bool isMatch;

	std::vector<double> tempAmpTimeVec;
	std::vector<double> tempAmpDataVec;

	std::smatch m;
	std::regex e("", std::regex_constants::icase);
	std::string s = "";
	gsl_matrix * ampArray;
	int nodeCnt = 0;

	std::istringstream ampStream(ampStringIn);
	while (std::getline(ampStream, ampLine, '\n'))
	{
		cnt++;
		if (cnt == 1)
		{
			//Get amplitude name
			e = std::regex("name=(.*?),", std::regex_constants::icase);
			s = ampLine;
			isMatch = std::regex_search(s, m, e);
			if (isMatch) sAmpName = m[1];
			else
			{
				e = std::regex("name=(.*)", std::regex_constants::icase);
				s = ampLine;
				isMatch = std::regex_search(s, m, e);
				if (isMatch) sAmpName = m[1];
				else
				{
					std::cout << "amplitude creation failed";
					return false;
				}
			}
		}
		else
		{
			int cnt1 = 0;
			std::istringstream ampDataStream(ampLine);
			while (std::getline(ampDataStream, ampDataTemp, ','))
			{
				if (cnt1 % 2 == 0) tempAmpTimeVec.push_back(std::stod(ampDataTemp));
				else tempAmpDataVec.push_back(std::stod(ampDataTemp));
				cnt1++;
			}
		}
	}

	if (tempAmpDataVec.size() < 2)
	{
		std::cout << "amplitude creation failed";
		return false;
	}

	if (tempAmpDataVec.size() != tempAmpTimeVec.size())
	{
		std::cout << "amplitude creation failed";
		return false;
	}

	ampArray = gsl_matrix_alloc(tempAmpDataVec.size(),2);
	std::vector<int>::const_iterator iter;
	cnt = 0;
	
	for (int i=0;i< tempAmpDataVec.size();i++)
	{
		gsl_matrix_set(ampArray, i, 0, tempAmpTimeVec[i]);
		gsl_matrix_set(ampArray, i, 1, tempAmpDataVec[i]);
	}	
	amplitudeMap[sAmpName] = ampArray;
	return true;
			

	//return false;
}

ModelData::~ModelData()
{
}

HistoryData::HistoryData()
{
}

bool HistoryData::addSteps(std::string szStepIn)
{
	std::string stepLine;
	std::string stepDataTemp;
	
	//std:string nodeset;
	int cnt = 0;
	std::string stepName,stepType;
	bool isMatch;	
	std::vector<double> tempDataVec;
	std::string stepBlockStr;
	std::smatch m;
	std::regex e("(^\\*step.*\\n[\\s\\S]*?)(?=\\n\\*)", std::regex_constants::icase);
	std::string s = szStepIn;
	isMatch = std::regex_search(s, m, e);
	if (isMatch) stepBlockStr = m[1];
	else
	{
		std::cout << "step creation failed";
		return false;
	}
	
	//----------------
	//Retrieve the step block data
	//----------------
	std::istringstream stepLineStream(stepBlockStr);
	while (std::getline(stepLineStream, stepLine, '\n'))
	{
		cnt++;
		if (cnt == 1)
		{
			//Get step name
			e = std::regex("name=(.*?),", std::regex_constants::icase);
			s = stepLine;
			isMatch = std::regex_search(s, m, e);
			if (isMatch) stepName = m[1];
			else
			{
				e = std::regex("name=(.*)", std::regex_constants::icase);
				s = stepLine;
				isMatch = std::regex_search(s, m, e);
				if (isMatch) stepName = m[1];
				else
				{
					std::cout << "step creation failed";
					return false;
				}
			}
		}
		else if (cnt == 2)
		{
			e = std::regex("^\\*(.*)", std::regex_constants::icase);
			s = stepLine;
			isMatch = std::regex_search(s, m, e);
			if (isMatch) stepType = m[1];
			else
			{
				std::cout << "step creation failed";
				return false;
			}
		}
		else // all time parameters
		{			
			std::istringstream stepDataStream(stepLine);
			while (std::getline(stepDataStream, stepDataTemp, ','))
			{
				tempDataVec.push_back(std::stod(stepDataTemp));
			}
		}
	}

	if (tempDataVec.size() == 0)
	{
		std::cout << "step creation failed";
		return false;
	}
		
	Step currStep(stepName);
	currStep.stepType = stepType;
	currStep.timeParameters = gsl_vector_alloc(tempDataVec.size());
	for (int i = 0; i < tempDataVec.size(); i++)
	{
		gsl_vector_set(currStep.timeParameters, i, tempDataVec[i]);
	}
	


	//----------------
	//Retrieve the BCs 
	//----------------
	s = szStepIn;
	try {
		std::regex elsetRegex("(^\\*boundary.*\\n[\\s\\S]*?)(?=\\n\\*)", std::regex_constants::icase);
		std::sregex_iterator next(s.begin(), s.end(), elsetRegex);
		std::sregex_iterator end;
		while (next != end) {
			m = *next;
			std::cout << m.str() << "\n";
			currStep.addBCtoSequence(m[1]);
			next++;
		}
	}
	catch (std::regex_error& err) {
		// Syntax error in the regular expression
		std::cout << "No node sets detected";
	}

	//----------------
	//Retrieve the Loads 
	//----------------
	s = szStepIn;
	try {
		std::regex elsetRegex("(^\\*cload.*\\n[\\s\\S]*?)(?=\\n\\*)", std::regex_constants::icase);
		std::sregex_iterator next(s.begin(), s.end(), elsetRegex);
		std::sregex_iterator end;
		while (next != end) {
			m = *next;
			std::cout << m.str() << "\n";
			currStep.addLoadtoSequence(m[1]);
			next++;
		}
	}
	catch (std::regex_error& err) {
		// Syntax error in the regular expression
		std::cout << "No node sets detected";
	}

	stepSequence.push_back(currStep);



	return true;
}

HistoryData::~HistoryData()
{
}

Step::Step(std::string stepNameIn)
{
	stepName = stepNameIn;
}

bool Step::addBCtoSequence(std::string bcStringIn)
{

	std::string bcLine;
	std::string bcDataLine;
	std::string bcDataTemp;
	//std:string nodeset;
	int cnt = 0;
	std::string bcName, bcType;
	bool isMatch;
	std::vector<std::string> tempDataVec;	
	std::smatch m;
	std::regex e("", std::regex_constants::icase);
	std::string s = "";	
	//bool isAmplitude =false;
	//std::string amplitudeName;

	BC currBC;
	currBC.isAmplitudeON = false;
	currBC.amplitudeName = "";
	currBC.BCType = "boundary";
	
	std::istringstream bcLineStream(bcStringIn);
	while (std::getline(bcLineStream, bcLine, '\n'))
	{
		cnt++;
		if (cnt == 1)
		{
			//Get amplitude
			e = std::regex("amplitude=(.*?),", std::regex_constants::icase);
			s = bcLine;
			isMatch = std::regex_search(s, m, e);
			if (isMatch) currBC.amplitudeName = m[1];
			else
			{
				e = std::regex("amplitude=(.*)", std::regex_constants::icase);
				s = bcLine;
				isMatch = std::regex_search(s, m, e);
				if (isMatch)
				{
					currBC.isAmplitudeON = true;
					currBC.amplitudeName = m[1];

				}				
			}
		}
		
		else // data lines
		{
			currBC.addtobcData(bcLine);
		}
	}

	BCs.push_back(currBC);

	return true;
}

bool Step::addLoadtoSequence(std::string loadStringIn)
{
	std::string loadLine;
	std::string loadDataLine;
	std::string loadDataTemp;

	//std:string nodeset;
	int cnt = 0;
	std::string loadName, loadType;
	bool isMatch;
	std::vector<std::string> tempDataVec;
	std::smatch m;
	std::regex e("", std::regex_constants::icase);
	std::string s = "";
	//bool isAmplitude =false;
	//std::string amplitudeName;

	Load currLoad;
	currLoad.isAmplitudeON = false;
	currLoad.amplitudeName = "";
	currLoad.loadType = "cload";

	std::istringstream loadLineStream(loadStringIn);
	while (std::getline(loadLineStream, loadLine, '\n'))
	{
		cnt++;
		if (cnt == 1)
		{
			//Get amplitude
			e = std::regex("amplitude=(.*?),", std::regex_constants::icase);
			s = loadLine;
			isMatch = std::regex_search(s, m, e);
			if (isMatch) currLoad.amplitudeName = m[1];
			else
			{
				e = std::regex("amplitude=(.*)", std::regex_constants::icase);
				s = loadLine;
				isMatch = std::regex_search(s, m, e);
				if (isMatch)
				{
					currLoad.isAmplitudeON = true;
					currLoad.amplitudeName = m[1];

				}
			}
		}

		else // data lines
		{
			currLoad.addLoadDataToSequence(loadLine);
		}
	}

	loads.push_back(currLoad);

	return true;

	
	
}

Step::~Step()
{
}


Load::Load()
{
}

bool Load::addLoadDataToSequence(std::string loadDataLineIn)
{
	LoadData currLoadData;
	std::vector<std::string> temp;
	std::string loadDataTemp;
	std::istringstream loadLineStream(loadDataLineIn);
	int cnt = 0;
	while (std::getline(loadLineStream, loadDataTemp, ','))
	{
		temp.push_back(loadDataTemp);
		cnt++;
	}
	if (cnt<3) return false;

	currLoadData.isNodeSetUsed = false;

	std::string t0;
	t0 = temp[0];


	//Check whether node set is used
	for (size_t j = 0; j < t0.length(); j++)
	{
		if (!isdigit(t0[j]))
			currLoadData.isNodeSetUsed = true;
		break;
	}

	if (currLoadData.isNodeSetUsed) currLoadData.nodeSet = t0;
	else currLoadData.nodeId = std::stoi(t0);

	currLoadData.dof = std::stoi(temp[1]);
	currLoadData.mag = std::stod(temp[2]);
	
	loadDataVec.push_back(currLoadData);
	return true;
}

Load::~Load()
{
}

BC::BC()
{
}

bool BC::addtobcData(std::string bcDataLineIn)
{
	BCData currBCData;
	std::vector<std::string> temp;
	std::string bcDataTemp;
	std::istringstream bcLineStream(bcDataLineIn);
	int cnt = 0;
	while (std::getline(bcLineStream, bcDataTemp, ','))
	{
		temp.push_back(bcDataTemp);
		cnt++;
	}
	if (cnt<3) return false;
		
	currBCData.isNodeSetUsed = false;

	std::string t0;
	t0 = temp[0];


	//Check whether node set is used
	for (size_t j = 0; j < t0.length(); j++)
	{
		if (!isdigit(t0[j]))
			currBCData.isNodeSetUsed = true;
			break;
	}

	if (currBCData.isNodeSetUsed) currBCData.nodeSet = t0;
	else currBCData.nodeId = std::stoi(t0);

	currBCData.startDOF = std::stoi(temp[1]);
	currBCData.endDOF = std::stoi(temp[2]);
	if (cnt == 4) currBCData.mag = std::stod(temp[3]);	
	bcData.push_back(currBCData);
	return true;
}

BC::~BC()
{
}

ResultData::ResultData()
{
}

bool ResultData::setSVec(int vecSize)
{
	std::vector<double> sVec; ;

	for (int i = 0; i < vecSize; i++)
	{
		sVec.push_back(0.0);
	}

	fields["S"] = sVec;
	return true;
}

ResultData::~ResultData()
{
}
