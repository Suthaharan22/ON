#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <vector>
#include <queue>
#include <set>
#include <algorithm>
#include <limits>
#include <assert.h>
#include <assert.h>
#include <time.h>

//=======================CHANGES WHEN GRAPH CHANGES================================
#define MAX_SLOTS 88	//Total number of slots 352
#define FIXED 50		//For a single wavelength (in GHz)
#define max_edges 42	//Changes when graph changes
#define min_nodes 1			//Minimum nodes
#define max_nodes 14			//Maximum nodes

#define TOTAL_PATHS 5		//Total paths we want to get (in allpathsmap) Always>=MAX_TRIAL_PATHS
#define MAX_TRIAL_PATHS 5			//Changes when graph changes (Total posible paths for primary and backup)
//=======================CHANGES WHEN GRAPH CHANGES================================

//=======================CHANGES WHEN LAYERS CHANGE================================
#define MAX_LAYERS 30			//Maximum size of layers
#define LAYERID 1			//Start value of Layer id
//=======================CHANGES WHEN LAYERS CHANGE================================

#define SHLPIDSTART 10001		//Initial Shared ID of a shared_id map

//=========================SIMULATION SETUP=====================================
#define MAX_LIGHTPATHS 100	//Total number of Lightpaths
#define DAT_RAT 100
//#define MAX_TRAFFIC	20		//Maximum traffic accomodating a network

#define MeanHoldTime 1		//IHT is Inter Holding Time (Mean). Usually it is 1.

#define LOAD10 10		//Load. and it varies.
#define LOAD20 20		//Load. and it varies.
#define LOAD50 500		//Load. and it varies.

#define TOT_LOAD 999999	//Total load for calculate sec eff
//=========================SIMULATION SETUP=====================================

using namespace std;

typedef vector<int> IntVec;
typedef map<int, IntVec> VecMap;
typedef multimap<int, IntVec> VecMul;
typedef set<int> IntSet;
typedef map<int, set<int>> MapSet;
typedef multimap<int, MapSet> TriMul;

//=========================FOR TRACKING CURRENT TIME STARTS===========================
time_t tim;
//=========================FOR TRACKING CURRENT TIME ENDS===========================

//=========================FUNCTION PROTOTYPES START===========================
//**************************FOR COPYING STARTS**************************
void insert_links();
void copy_slots_mylink_mainmap();
void display_mainmap();
//**************************FOR COPYING ENDS**************************

//**************************FOR GRAPH STARTS**************************
void read_graph_from_text();
void display_graph_text();
void display_fiberid_wt_text();
void display_fiber_link_text();

void MoreThanFourLinks(VecMap& graphtextmap);
void display_MoreThanFourLinks();
//**************************FOR GRAPH ENDS**************************

//**************************FOR DATA RATE STARTS**************************
void get_data_rate();
void display_data_rate();
int read_data_rate(int linerate);
//**************************FOR DATA RATE STARTS**************************

//-------------------------FOR FIND ALL PATHS STARTS---------------------
bool isadjacency_node_not_present_in_current_path(int node, IntVec pathway);
int findallpaths(int source, int target);
void display_all_paths_mmap();

void arrangement_of_links();
void display_arranged_links();

void collection_of_paths();
void display_pathmaps();

//-------------------------FOR FIND ALL PATHS STARTS---------------------

//-------------------------START CONVERTING NODES TO LINKS---------------------------
vector<int> convert_to_links(const IntVec& nodevec);
//-------------------------END CONVERTING NODES TO LINKS---------------------------

//------------------------FOR FIND BACKUP PATH (1) STARTS-------------------------
bool imatch(const IntVec& vec, int x, int y);
bool vmatch(const IntVec& vec1, const IntVec& vec2);
//------------------------FOR FIND BACKUP PATH (1) ENDS-------------------------

//*************************FOR SELECTING PRIMARY***************************
void NEW_select_primary_path();		//For selecting primary slots and selecting backup paths
bool PrimaryAllAllowRange(const IntVec& pplinks, int startIndex, int len, const VecMap& theMap);
bool Allow(int value);
//*************************FOR SELECTING PRIMARY***************************

//------------------------FOR FIND BACKUP PATH (2) STARTS-------------------------
bool Find(const IntVec& vec, int value);
bool HaveNoCommonValueEither(const map<int, int>& themap, const IntVec& vec1, const IntVec& vec2);
bool AllowKey(const VecMap& primap, int specialKey, int key);
VecMap CreateNewMap(const VecMap& primap, int specialKey, const VecMap& secmap, const IntVec& specialVec);
void display_sorted_map();
//------------------------FOR FIND BACKUP PATH (2) ENDS-------------------------

//*************************FOR INSERTING SPECIAL SHAREABLE IDS**********************************
bool IsValueInMap(int value, const VecMap& m2);
bool AreAllValuesInMap(const IntVec& v, const VecMap& m2);
IntVec other_shareable_ids(const VecMap& sorted_map);
//************************FOR INSERTING SPECIAL SHAREABLE IDS**********************************

//*************************FOR SELECTING SECONDARY PATH STARTS***************************
bool SecondaryAllAllowRange(const IntVec& mvec, int startIndex, int len, const VecMap& theMap, const VecMap& m1, const IntVec& v);
bool SecondaryAllow(int value, int key, const VecMap& m1, const IntVec& v);
int SecondaryCountAllSpecialInRange(const IntVec& mvec, int startIndex, int len, VecMap& theMap, const VecMap& m1, const IntVec& v);
bool SecondarySpecial(int value, const VecMap& m1, const IntVec& v);

VecMul CreateNewVector(const IntVec& myvec, VecMap& mainmap, const VecMap& newmap, const IntVec& vect);
void display_slots_range();			//For Displaying slots range with their counts
//*************************FOR SELECTING SECONDARY PATH ENDS***************************

//*************************FOR INSERTING WHOLE PRIMARY & SECONDARY PATH STARTS***************************
int FindMapKeyWithSameValues(VecMap& smap, IntVec& vec);
void ProcessData(IntVec& mVec, VecMap& smap, int j, int positionKey);

void select_slots_for_secondarypath();

void insert_all_lightpaths();
//*************************FOR INSERTING WHOLE PRIMARY & SECONDARY PATH ENDS***************************

//*************************FOR INSERTING PRIMARY AND BACKUP PATHS*****************************************
void insert_backup_lp(const IntVec& backvec);
void insert_primary_lp(IntVec& primvec);
//*************************FOR INSERTING PRIMARY AND BACKUP PATHS*****************************************

//**************************DISPLAYING STARTS**************************
void display_pplinks();

void display_pslotposition();
void display_primarypath();
void display_primarynode();

void display_bslotposition();
void display_backuppath();
void display_backupnode();

void display_shared_id();	//display shared_id and lightpaths

void display_deleted_lightpaths();	//displaying deleted lightpaths from a vector
//**************************DISPLAYING ENDS**************************

//*************************FOR TRIAL STARTS***************************
void admit_lightpath(int start, int end, int Lightpath);
void depart_lightpath(int dLightpath);

bool lpchecker(VecMap& mypbmaps);
void printcase1();
void printcase2();
void printcase3();
void printcase4();
void printcase5();
void printcase6();
//*************************FOR TRIAL ENDS***************************

//*************************FOR DELETING PRIMARY & SECONDARY PATH STARTS***************************
void del_pri_lp_from_mainmap(int dLightpath);
void ProcessDeleteData(IntVec& mVec, int j, int positionKey);

void del_back_lp_from_mainmap(int dLightpath);
bool Find(const vector<int>& myvector, const map<int, vector<int>>& vecMap, int value, const pair<int, int>& position);
int NewValue(int key, const vector<int>& vec, const map<int, vector<int>>& smap);
void Replace(const vector<int>& dellinks, map<int, vector<int>>& mainmap, int oldValue, int newValue, const pair<int, int>& position);
//*************************FOR DELETING PRIMARY & SECONDARY PATH ENDS***************************

//=========================IMPLEMENTING FLEX-WSS STARTS=================================
//**************************FOR COPYING (FOR FLEX-WSS) STARTS**************************
void insert_nodes();
void copy_layers_mynode_mainnodemap();
void display_mainnodemap();
//**************************FOR COPYING (FOR FLEX-WSS) ENDS**************************

//-----------------------MAIN FLEX-WSS IMPLEMENTATION-------------------------------------
map<int, int> ProcessLayer(int Sslot, int Eslot);
//-----------------------MAIN FLEX-WSS IMPLEMENTATION-------------------------------------

//---------------------FOR VARIOUS CASE LINK OPERATION------------------------------------
map<int, int> BothNullLinkOperation(int node);
map<int, int> AnyOneLinkOperation(IntVec& AnyOneLink, int node);
map<int, int> BothFullLinkOperation(IntVec& BeforeVec, IntVec& AfterVec, int node);
//---------------------FOR VARIOUS CASE LINK OPERATION------------------------------------

//------------------------SUPPORTING FUNCTIONS----------------------------------------------
bool IsValueInVec(int value, IntVec& v2);
bool isVectorMatch(IntVec& v1, IntVec& v2);
bool check_current_lp_overlap(int Sslot, int Eslot, int lp);
//------------------------SUPPORTING FUNCTIONS----------------------------------------------

//------------------------FOR GENERATING POSSIBLE SHARABLE LINKS-------------------------------------
IntVec BeforeAfterLinks(VecMap& sorted_map, int link, int Sslot, int Eslot);
//------------------------FOR GENERATING POSSIBLE SHARABLE LINKS-------------------------------------

//------------------FOR DISPLAYING NODE_LAYER FOR CURRENT LP-----------------------------------------
void display_node_layer();	//Displaying stored node and layer map	 for current lp
//------------------FOR DISPLAYING NODE_LAYER FOR CURRENT LP-----------------------------------------

//------------------FOR INSERTING DELETING DISPLAYING LIGHTPATH INTO MAINNODEMAP AND LAYER_ID------------
void insert_lightpath_into_layer();
void del_lightpath_in_layers(int dLightpath);
void displaying_layer_id();
//------------------FOR INSERTING DELETING DISPLAYING LIGHTPATH INTO MAINNODEMAP AND LAYER_ID------------

//------------------FOR INSERTING DELETING DISPLAYING LIGHTPATH IN NODE_LPS-----------------------------
void insert_lightpath_into_node_lps();
void del_lightpath_in_node_lps(int dLightpath);
void display_node_lps();
//------------------FOR INSERTING DELETING DISPLAYING LIGHTPATH IN NODE_LPS-----------------------------
//=========================IMPLEMENTING FLEX-WSS ENDS=================================

//********************************FOR PBPS STARTS**************************************************
//------------------FOR INSERTING, DISPLAYING, DELETING LPS IN BPLINK_LPS-----------------------------
void insert_lightpath_into_link_lps();
void displaying_bplink_lps();
void del_lightpath_from_bplink_lps(IntVec& bplVec, int dLightpath);
//void del_lightpath_from_bplink_lps_1(int Lightpath);
//------------------FOR INSERTING, DISPLAYING, DELETING LPS FROM BPLINK_LPS-----------------------------

//------------------FOR INSERTING, DISPLAYING, DELETING LPS IN SPLITMULMAP-----------------------------
void insert_splitmulmap(const IntVec& backvec);
bool NotFind(set<int>& myset, int value);
void display_splitmulmap();
void del_lightpath_from_splitmulmap(int dLightpath);
//------------------FOR INSERTING, DISPLAYING, DELETING LPS IN SPLITMULMAP-----------------------------

//==============================For Checking current lp with other history (CASE-1)===========================================
bool check_is_curr_joint(IntVec& curVec, int lp);
bool check_current_lp_overlap(int Sslot, int Eslot, int lp);
//==============================For Checking current lp with other history (CASE-1)========================================
//===========================For Checking other lps with current lp history (CASE-2)==========================================
bool check_is_joint(int lp1, int lp2);
bool check_overlap_bp_slots(int lp1, int lp2);
//============================For Checking other lps with current lp history (CASE-2)================================

//=============================MAIN PBPS FUNCTION====================================
IntSet PBPS_COND_1(int start, int end);
bool FindSet(const IntSet& iSet, int value);
IntSet PBPS_COND_2(int start, int end);
//=============================MAIN PBPS FUNCTION====================================
//********************************FOR PBPS ENDS**************************************************

//=============================FOR AUTOMATION STARTS====================================
int random_source(int lowest_number, int highest_number);	//Creating source node
int random_destination(int lowest_number, int highest_number, int Source);	//Creating destination node
double expon(double mean);			//returns random exponential-variate with specified mean
void display_TIME_SCALE();
//=============================FOR AUTOMATION ENDS====================================

//===========================RESULTS==========================================
void Calculation();
//===========================RESULTS==========================================
//=======================FUNCTION PROTOTYPES END===========================


//==================================DECLARATIONS===========================================
//--------------------------------SETTING LINKS START--------------------------------------------------
map<int, int*> mylink;				//For inserting links and slots temporarily
map<int, int*>::iterator itl;

VecMap mainmap;				//For inserting links and slots permanently
VecMap::iterator itmm;
//---------------------------------SETTING LINKS END-------------------------------------------------

//----------------------------FOR GRAPH STARTS------------------------------------------------
VecMap graphtextmap;			//For inserting path and it's Link_ID from Link_ID text file
VecMap::iterator itgtm;

IntSet FourLinks;				//For storing links in which that link is one of the 4 or more connection with the node.

IntVec pplinks;					//pplinks - Primary Path Links
IntVec::iterator itppl;

IntVec ppnodes;					//ppnodes - Primary Path Nodes
IntVec::iterator itppn;

IntVec bplinks;					//bplinks - Backup Path Links
IntVec::iterator itbpl;

IntVec bpnodes;					//bpnodes - Backup Path Nodes
IntVec::iterator itbpn;

map<int, int> fiberid_weight_map;
map<int, int>::iterator fbwtmp;

vector<vector<int> >GRAPH(100);

IntVec pathway;

VecMul allpathsmap;			//For inserting all possible paths from source to destination
VecMul::iterator itapm;

VecMul pathsmmap;			//For inserting possible links with their weight
VecMul::iterator pmmap;

VecMul arranged_mmap;			//For inserting possible links with their weight
VecMul::iterator armmap;
//----------------------------FOR GRAPH ENDS------------------------------------------------

//---------------------------FOR FIBER<=>LINKS STARTS-------------------------------------------
map<int, int> fiber_link;			//For inserting fiber vs. links
map<int, int>::iterator itfl;
//---------------------------FOR FIBER<=>LINKS STARTS---------------------------------------------

//---------------------------FOR DATA RATE STARTS-------------------------------------------
multimap<int, int> DataRatemm;			//For inserting data rate
multimap<int, int>::iterator itdrmm;
//---------------------------FOR DATA RATE ENDS---------------------------------------------

//----------------------------FOR INSERTING LIGHTPATHS START------------------------------------------
VecMap primarypath;				//For inserting paths and links of primary paths
VecMap::iterator itppath;

VecMap primarynode;				//For inserting paths and nodes of primary paths
VecMap::iterator itpnode;

VecMap backuppath;				//For inserting paths and links of backup paths
VecMap::iterator itbpath;

VecMap backupnode;				//For inserting paths and nodes of backup paths
VecMap::iterator itbnode;
//----------------------------FOR INSERTING LIGHTPATHS END--------------------------------------------

//---------------------------FOR SLOTS STARTS-------------------------------------------
VecMap prislotposition;				//For inserting primary paths and their corresponding slot positions
VecMap::iterator psltposn;

VecMap backslotposition;				//For inserting backup paths and their corresponding slot positions
VecMap::iterator bsltposn;
//---------------------------FOR SLOTS ENDS---------------------------------------------

//---------------------------FOR SHARING STARTS-------------------------------------------
VecMap sorted_map;		//For inserting sharable lightpaths

VecMap shared_id;		//For inserting shared lightpaths ID and lightpaths
VecMap::iterator slpsid;

IntVec splvec;			//To insert special values from shared_id comparing sorted_map

VecMul VMap;		//To insert selected range of slots with their count
//---------------------------FOR SHARING ENDS--------------------------------------------

//---------------------------FOR FLEX-WSS STARTS-------------------------------------------
map<int, int*> mynode;				//For inserting nodess and layes temporarily
map<int, int*>::iterator itn;

VecMap mainnodemap;				//For inserting nodes and layers permanently
VecMap::iterator mnodemap;

VecMap layer_id;				//For inserting layer_id and lightpaths

map<int, set<int>> node_lps;		//For storing nodes and sharable lightpaths for the current lightpath

map<int, int> NodeLayer;		//For inserting nodes and layers for a particular lightpath
//---------------------------FOR FLEX-WSS ENDS-------------------------------------------

//--------------------------FOR PBPS STARTS---------------------------------------------
map<int, set<int>> bplink_lps;
TriMul SplitMulMap;

IntSet NotifySet1;
IntSet NotifySet2;
//--------------------------FOR PBPS ENDS---------------------------------------------

//--------------------------FOR DELETE FUNCTIONS----------------------------------------
IntVec DelVec;		//Inserting deleted lightpaths
IntVec::iterator dVec;
//--------------------------FOR DELETE FUNCTIONS----------------------------------------

//--------------------------FOR DISTRIBUTION----------------------------------------
multimap<double, int> TIME_SCALE;		//Inserting time scale of the lightpaths
//--------------------------FOR DISTRIBUTION----------------------------------------

//--------------------------FOR COUNTING BLOCKING------------------------------------
IntSet NOPRIPATH;
IntSet NOPRISLOT;
IntSet NOBACKPATH;
IntSet NOBACKSLOT;
IntSet NOLAYER;
IntSet NOADDCON;
IntSet XBLOCK;

IntVec NOBACKSLOTV;
IntVec NOLAYERV;
IntVec NOADDCONV;
IntVec XBLOCKV;

IntVec RATETOTAL;
IntVec RATESBLOCKED;
//--------------------------FOR COUNTING BLOCKING------------------------------------

int start, stop;
int Lightpath, Data_Rate, No_Of_Slots;
int dLightpath;

double MeanArrivalTime, aTime, hTime;
double curTime = 0;

int LOAD;

ofstream fout("PBPS_Results.txt", ios::out | ios::app);

//STARTING OF MAIN FUNCTION
int main(int argc, const char* argv[])
{

	//=========================FOR TRACKING CURRENT TIME STARTS===========================
	time(&tim);
	char mytime[123];
	ctime_s(mytime, _countof(mytime), &tim);
	cout << "\n----------------------START OF SIMULATION---------------------\n";
	cout << "\nSimulation Start Time: " << mytime;
	fout << "\n----------------------START OF SIMULATION---------------------\n";
	fout << "\nSimulation Start Time: " << mytime;
	//=========================FOR TRACKING CURRENT TIME STARTS===========================

	//-------------------------INSERTING SLOTS INTO MAINMAP----------------------------------
	insert_links();		// insert total slots to the total edges
	copy_slots_mylink_mainmap();	//copy all the links' slots to the mainmap

		//display_mainmap();

	//-------------------------INSERTING SLOTS INTO MAINMAP----------------------------------

	//-------------------------INSERTING LAYERS INTO MAINNODEMAP----------------------------------
	insert_nodes();
	copy_layers_mynode_mainnodemap();
	//display_mainnodemap();		//display the whole nodes and their layers which are in mainnodemap map
//-------------------------INSERTING LAYERS INTO MAINNODEMAP----------------------------------

	read_graph_from_text();		//create a graph
		//display_graph_text();	//display graph from Link_ID.txt->multimap->display
		//display_fiberid_wt_text();
		//display_fiber_link_text();

	MoreThanFourLinks(graphtextmap);	//For storing Link which node connected with >4 links
		//display_MoreThanFourLinks();	//Displaying Link which node connected with >4 links

	get_data_rate();
	//display_data_rate ();															


//**************************TO RUN MULTIPLE LOADS**********************************************
	IntVec LOADVec;

	//LOADVec.push_back(LOAD10); 
	//LOADVec.push_back(LOAD20); 
	LOADVec.push_back(LOAD50);

	for (IntVec::iterator itload = LOADVec.begin(); itload != LOADVec.end(); ++itload)
	{
		LOAD = (*itload);

		//---------------------INITIALIZING ALL CONTAINERS START---------------------------------------
		//mainmap.clear();				//For inserting links and slots permanently
		for (VecMap::iterator itvmap = mainmap.begin(); itvmap != mainmap.end(); ++itvmap)
		{
			IntVec& vecmap = itvmap->second;
			for (size_t j = 0; j < vecmap.size(); j++)
			{
				if (vecmap[j] != 0)
				{
					vecmap[j] = 0;
				}
			}
		}
		//----------------------------------------------------------------------------------------------
		//------------------------------------------------------------------------------------------------
		primarypath.clear();				//For inserting paths and links of primary paths
		primarynode.clear();				//For inserting paths and nodes of primary paths
		backuppath.clear();				//For inserting paths and links of backup paths
		backupnode.clear();				//For inserting paths and nodes of backup paths
		//-----------------------------------------------------------------------------------------------
		//------------------------------------------------------------------------------------------------
		prislotposition.clear();				//For inserting primary paths and their corresponding slot positions
		backslotposition.clear();				//For inserting backup paths and their corresponding slot positions
		//------------------------------------------------------------------------------------------------
		//--------------------------------------------------------------------------------------------------
		shared_id.clear();		//For inserting shared lightpaths ID and lightpaths
		//----------------------------------------------------------------------------------------------------
		//----------------------------------------------------------------------------------------------------
		//mainnodemap.clear();				//For inserting nodes and layers permanently
		for (VecMap::iterator itnodemap = mainnodemap.begin(); itnodemap != mainnodemap.end(); ++itnodemap)
		{
			IntVec& nodemap = itnodemap->second;
			for (size_t k = 0; k < nodemap.size(); k++)
			{
				if (nodemap[k] != 0)
				{
					nodemap[k] = 0;
				}
			}
		}

		layer_id.clear();				//For inserting layer_id and lightpaths
		node_lps.clear();		//For storing nodes and sharable lightpaths for the current lightpath
		//---------------------------------------------------------------------------------------------------
		//----------------------------------------------------------------------------------------------------
		bplink_lps.clear();
		SplitMulMap.clear();
		//-----------------------------------------------------------------------------------------------------
		//----------------------------------------------------------------------------------------------------
		DelVec.clear();		//Inserting deleted lightpaths
		//----------------------------------------------------------------------------------------------------
		//----------------------------------------------------------------------------------------------------
		TIME_SCALE.clear();		//Inserting time scale of the lightpaths
		//----------------------------------------------------------------------------------------------------
		//----------------------------------------------------------------------------------------------------
		NOPRIPATH.clear();
		NOPRISLOT.clear();
		NOBACKPATH.clear();
		NOBACKSLOT.clear();
		NOLAYER.clear();
		NOADDCON.clear();
		XBLOCK.clear();

		RATETOTAL.clear();
		RATESBLOCKED.clear();
		//----------------------------------------------------------------------------------------------------
		//========================INITIALIZING ALL CONTAINERS ENDS=================================

		srand(static_cast<unsigned>(time(0)));

		MeanArrivalTime = (double)MeanHoldTime / LOAD;

		int Rate10 = 0;
		int Rate40 = 0;
		int Rate100 = 0;
		int Rate400 = 0;
		int Rate1000 = 0;

		for (int i = 1; i <= MAX_LIGHTPATHS; i++)
		{
			/*
					/*
					cout<<"Enter the Lightpath: ";		//<<Lightpath;
					cin>>Lightpath;

					cout<<"\nEnter the Data Rate to be used: ";		//<<No_Of_Slots;
					cin>>Data_Rate;

					read_data_rate();		//getting No_Of_Slots according to the Data_Rate

					cout<<"\nEnter source node: ";		//<<start;
					cin>>start;

					cout<<"\nEnter destination node: ";		//<<end;
					cin>>end;
			*/

			//--------------Automated----------------------------------
			//---------------For Lightpath----------------------------
			Lightpath = i;
			cout << "\nLightpath to be admitted: " << Lightpath;
			//---------------For Lightpath----------------------------

			//---------------For Data Rate----------------------------
			multimap<int, int>::iterator itmul = DataRatemm.begin();
			advance(itmul, rand() % DataRatemm.size());
			int random_key = itmul->first;

			Data_Rate = random_key;
			cout << "\nData Rate: " << Data_Rate << " Gbps";

			//For counting total rates of total lightpaths
			if (Data_Rate == 10)
			{
				Rate10++;
			}
			else if (Data_Rate == 40)
			{
				Rate40++;
			}
			else if (Data_Rate == 100)
			{
				Rate100++;
			}
			else if (Data_Rate == 400)
			{
				Rate400++;
			}
			else if (Data_Rate == 1000)
			{
				Rate1000++;
			}

			No_Of_Slots = read_data_rate(Data_Rate);		//getting No_Of_Slots according to the Data_Rate

			cout << "\nFrequency Slots for Lightpath: " << No_Of_Slots;
			//---------------For Data Rate----------------------------

			//---------------For Source Node----------------------------
			start = random_source(min_nodes, max_nodes);
			cout << "\nSource Network Node: " << start;
			//---------------For Source Node----------------------------

			//---------------For Destination Node----------------------------
			stop = random_destination(min_nodes, max_nodes, start);
			cout << "\nDestination Network Node: " << stop;
			//---------------For Destination Node----------------------------

			//---------------For Arrival Time----------------------------
			aTime = expon(MeanArrivalTime);
			cout << "\nArrival Time of a Lightpath: " << aTime;
			//---------------For Arrival Time----------------------------

			//---------------For Holding Time----------------------------
			hTime = expon(MeanHoldTime);
			cout << "\nHolding Time of a Lightpath: " << hTime;
			//---------------For Holding Time----------------------------
			//--------------Automated----------------------------------

			cout << "\n\n";

			//++++++++++++++++++++++++--------AUTOMATION----------++++++++++++++++++++++++++++++++++
			curTime += aTime;

			if (TIME_SCALE.empty())
			{
				admit_lightpath(start, stop, Lightpath);

				if (lpchecker(backslotposition))
				{
					cout << "\nCURRENT TIME: " << curTime;
					TIME_SCALE.insert(make_pair((curTime + hTime), Lightpath));
				}
			}
			else if (!TIME_SCALE.empty())
			{
				if (curTime <= TIME_SCALE.begin()->first)
				{
					admit_lightpath(start, stop, Lightpath);
					if (lpchecker(backslotposition))
					{
						cout << "\nCURRENT TIME: " << curTime;
						TIME_SCALE.insert(make_pair((curTime + hTime), Lightpath));
					}

				}
				else if (curTime > TIME_SCALE.begin()->first)
				{
					double delcurTime = 0;
					do
					{
						delcurTime = TIME_SCALE.begin()->first;
						cout << "\nCURRENT TIME: " << delcurTime;
						dLightpath = TIME_SCALE.begin()->second;

						multimap<double, int>::const_iterator itscale = TIME_SCALE.begin();
						TIME_SCALE.erase(itscale);

						depart_lightpath(dLightpath);

					} while ((!TIME_SCALE.empty()) && (curTime > delcurTime));

					admit_lightpath(start, stop, Lightpath);
					if (lpchecker(backslotposition))
					{
						TIME_SCALE.insert(make_pair((curTime + hTime), Lightpath));
						cout << "\nCURRENT TIME: " << curTime;
					}
				}
			}

			display_TIME_SCALE();		//Displaying time scale for the automation
	//++++++++++++++++++++++++--------AUTOMATION----------++++++++++++++++++++++++++++++++++

		}// End of for(int i=1; i<=MAX_LIGHTPATHS; i++)

		//Inserting total number of rates into a vector
		RATETOTAL.push_back(Rate10);
		RATETOTAL.push_back(Rate40);
		RATETOTAL.push_back(Rate100);
		RATETOTAL.push_back(Rate400);
		RATETOTAL.push_back(Rate1000);
		//Inserting total number of rates into a vector

		Calculation();	//Performance of the simulation and output to the text file

		cout << "\n";

	}	//for(IntVec::iterator itload=LOADVec.begin(); itload!=LOADVec.end();++itload)

	return 0;

}	//End of main function


//===================================INSERT MAINMAP STARTS=============================================
//INSERT LINKS AND SLOTS INITIALLY
void insert_links()
{
	for (int i = 1; i <= max_edges; i++)
	{
		mylink.insert(make_pair((i), new int[MAX_SLOTS]()));		//(*itslp) indicates the links of shortest path
	}
}

//COPYING SLOTS FROM MYLINK TO MAINMAP
void copy_slots_mylink_mainmap()
{
	for (itl = mylink.begin(); itl != mylink.end(); ++itl)
	{
		for (size_t i = 0; i < MAX_SLOTS; ++i)
		{
			mainmap[itl->first].push_back(itl->second[i]);
		}
	}
	mylink.erase(mylink.begin(), mylink.end());
}

//DISPLAY MAINMAP MULTIMAP
void display_mainmap()
{
	cout << "\n\n\nMAINMAP CONTAINS (LINK=>SLOTS)\n";
	cout << "------------------------------";
	for (itmm = mainmap.begin(); itmm != mainmap.end(); ++itmm)
	{
		cout << "\n" << (*itmm).first << "=>";
		for (size_t n = 0; n < (*itmm).second.size(); n++)
		{
			cout << (*itmm).second[n] << "  ";
		}
	}
	cout << "\n";
}
//=====================================INSERT MAINMAP ENDS=============================================

//=================================START OF DATA RATE FUNCTION========================================
//READ DATA RATE FROM TEXT FILE
void get_data_rate()
{
	int DR, NOS;
	ifstream file2("Data_Rate.txt");
	if (!file2.is_open())
	{
		cout << "Error opening file Data_Rate.txt"; return;
	}

	while (!file2.eof())
	{
		file2 >> DR >> NOS;
		DataRatemm.insert(make_pair(DR, NOS));
	}
	file2.close();

}

//DISPLAY DATA RATE FROM DATARATE MULTIMAP
void display_data_rate()
{
	for (itdrmm = DataRatemm.begin(); itdrmm != DataRatemm.end(); ++itdrmm)
	{
		cout << itdrmm->first << "=>" << itdrmm->second << "\n";
	}
}

//GETTING No_Of_Slots ACCORDING TO THE DATA RATE
int read_data_rate(int linerate)
{
	for (itdrmm = DataRatemm.begin(); itdrmm != DataRatemm.end(); ++itdrmm)
	{
		if (linerate == (*itdrmm).first)
		{
			return((*itdrmm).second);
		}
	}
	return 0;
}
//=================================END OF DATA RATE FUNCTION========================================

//==================================START OF GRAPH FUNCTIONS===========================================
//CREATE A GRAPH
void read_graph_from_text()
{
	int origin, destin, weight, fiber_id, link_id;

	//Link_ID
	ifstream file1("Link_ID.txt");
	if (!file1.is_open())
	{
		cout << "Error opening file Link_ID.txt"; return;
	}

	while (!file1.eof())
	{
		for (int index = 0; index < max_edges; index++)
		{
			file1 >> origin >> destin >> weight >> fiber_id >> link_id;

			GRAPH[origin].push_back(destin);				//FOR SURVIVABILITY....

			graphtextmap[fiber_id].push_back(origin);
			graphtextmap[fiber_id].push_back(destin);			//Inserting fiberid and source-destination into graphtextmap map from Link_ID text file

			fiberid_weight_map.insert(make_pair(fiber_id, weight));		//Inserting fiberid and weight into fiberid_weight_map map from Link_ID text file

			fiber_link.insert(make_pair(fiber_id, link_id));		//Inserting fiberid and linkid into fiber_link map from Link_ID text file	
		}
		cout << "\n";
	}

	file1.close();
}//End of create graph()

//DISPLAY GRAPH FROM TEXT FILE
void display_graph_text()
{
	cout << "\n\nGRAPHTEXT MAP (FIBERID=>SOURCE-DEST.)\n";
	cout << "------------------";
	for (itgtm = graphtextmap.begin(); itgtm != graphtextmap.end(); ++itgtm)
	{
		cout << endl << itgtm->first << " => ";
		for (size_t n = 0; n < (*itgtm).second.size(); ++n)
		{
			cout << " " << itgtm->second[n] << " ";
		}
	}
	cout << "\n";
}

//DISPLAY LINKID_WEIGHT FROM TEXT FILE
void display_fiberid_wt_text()
{
	cout << "\n\nFIBERID MAP (FIBERID=>WEIGHT)\n";
	cout << "------------------";
	for (fbwtmp = fiberid_weight_map.begin(); fbwtmp != fiberid_weight_map.end(); ++fbwtmp)
	{
		cout << "\n" << fbwtmp->first << "=>" << fbwtmp->second;
	}
}

//DISPLAY FIBER_LINK FROM TEXT FILE												
void display_fiber_link_text()
{
	cout << "\n\nFIBER_LINK MAP (FIBERID=>LINKID)\n";
	cout << "------------------";
	for (itfl = fiber_link.begin(); itfl != fiber_link.end(); ++itfl)
	{
		cout << "\n" << itfl->first << "=>" << itfl->second;
	}
}

//FOR STORING MORE THAN FOUR LINKS OF A NODE
void MoreThanFourLinks(VecMap& graphtextmap)
{
	VecMap tempM;
	tempM.clear();
	for (VecMap::iterator it = graphtextmap.begin(); it != graphtextmap.end(); ++it)
	{
		int key = it->first;
		IntVec& vec = it->second;
		for (size_t n = 0; n < vec.size(); ++n)
		{
			tempM[vec[n]].push_back(key);
		}
	}

	for (VecMap::iterator itm = tempM.begin(); itm != tempM.end(); ++itm)
	{
		IntVec& tempV = itm->second;
		if (tempV.size() >= 8)
		{
			for (IntVec::iterator itv = tempV.begin(); itv != tempV.end(); ++itv)
			{
				FourLinks.insert(*itv);
			}
		}
	}
}

//FOR DISPLAYING MORE THAN FOUR LINKS OF A NODE
void display_MoreThanFourLinks()
{
	cout << "\n\nMORE THAN FOUR LINKS\n";
	cout << "---------------------";
	for (IntSet::iterator it = FourLinks.begin(); it != FourLinks.end(); ++it)
	{
		cout << endl << (*it);
	}
	cout << "\n";
}
//**********************************END OF GRAPH FUNCTIONS*******************************************

//===================================FOR ALL PATHS STARTS===========================================
//ADJACENCY NODE NOT PRESENT IN CURRENT PATH
bool isadjacency_node_not_present_in_current_path(int node, vector<int>pathway)
{
	for (unsigned int i = 0; i < pathway.size(); ++i)
	{
		if (pathway[i] == node)
			return false;
	}
	return true;
}

//FINDS ALL PATHS		====================> allpathsmap
int findallpaths(int source, int target)
{
	int mykey = 0;

	pathway.clear();
	allpathsmap.clear();
	pathway.push_back(source);
	queue<vector<int> >q;
	q.push(pathway);

	while (!q.empty())
	{
		pathway = q.front();
		q.pop();

		int last_nodeof_path = pathway[pathway.size() - 1];
		if (last_nodeof_path == target)
		{

			int totwt = 0;
			for (size_t i = 1; i != pathway.size(); i++)
			{

				for (itgtm = graphtextmap.begin(); itgtm != graphtextmap.end(); ++itgtm)
				{
					const vector<int>& vgtm = itgtm->second;
					{
						if (((pathway[i - 1] == vgtm[0]) && (pathway[i] == vgtm[1])))
						{
							for (fbwtmp = fiberid_weight_map.begin(); fbwtmp != fiberid_weight_map.end(); ++fbwtmp)
							{
								if (itgtm->first == fbwtmp->first)
								{
									totwt = totwt + fbwtmp->second;

									//v.push_back(pathway[i-1]);
								}
							}
						}

					}
				}

			}
			allpathsmap.insert(make_pair(totwt, pathway));

			//-----FOR SELECTING ONLY CERTAIN NUMBER OF PATHS (NOT SHORTEST PATHS)-------
			if (allpathsmap.size() == TOTAL_PATHS)
			{
				break;
			}
			//-----FOR SELECTING ONLY CERTAIN NUMBER OF PATHS (NOT SHORTEST PATHS)-------	

		}
		else
			//print_path(path);

			for (unsigned int i = 0; i < GRAPH[last_nodeof_path].size(); ++i)
			{
				if (isadjacency_node_not_present_in_current_path(GRAPH[last_nodeof_path][i], pathway))
				{
					vector<int>new_path(pathway.begin(), pathway.end());
					new_path.push_back(GRAPH[last_nodeof_path][i]);
					q.push(new_path);
				}
			}
	}
	return 1;
}

//DISPLAY ALL POSSIBLE PATHS FROM MULTIMAP ALLPATHMM		
void display_all_paths_mmap()
{
	cout << "\n\nALL POSSIBLE PATHS\n";
	cout << "------------------";
	for (itapm = allpathsmap.begin(); itapm != allpathsmap.end(); ++itapm)
	{
		cout << endl << itapm->first << " => ";
		for (size_t n = 1; n < (*itapm).second.size(); ++n)
		{
			cout << itapm->second[n - 1] << "->" << itapm->second[n] << " ";
		}
	}
	cout << "\n";
}
//======================================FOR ALL PATHS ENDS========================================

//=========================SELECTING NUMBER OF PATHS FOR PRIMARY AND BACKUP=================================
//SELECT PATHS FOR TRIALS			====================> allpathsmap to pathsmmap
void collection_of_paths()
{
	pathsmmap.clear();
	vector<int> tempvec;
	for (itapm = allpathsmap.begin(); itapm != allpathsmap.end(); ++itapm)
	{
		if (pathsmmap.size() < MAX_TRIAL_PATHS)
		{
			for (size_t i = 0; i < (*itapm).second.size(); ++i)
			{
				tempvec.push_back(itapm->second[i]);
			}
			pathsmmap.insert(make_pair(itapm->first, tempvec));
			tempvec.clear();
		}
		if (pathsmmap.size() == MAX_TRIAL_PATHS)
		{
			break;
		}
	}
}

//DISPLAY PATHSMMAP FROM WEIGHTMAP MULTIMAP
void display_pathmaps()
{
	cout << "\n\nSELECTED PATHSMMAP\n";
	cout << "------------------";
	for (pmmap = pathsmmap.begin(); pmmap != pathsmmap.end(); ++pmmap)
	{
		cout << endl << pmmap->first << " => ";
		for (size_t n = 1; n < (*pmmap).second.size(); ++n)
		{
			cout << pmmap->second[n - 1] << "->" << pmmap->second[n] << " ";
		}
	}
}
//=========================SELECTING NUMBER OF PATHS FOR PRIMARY AND BACKUP=================================

//=====================ARRANGING LINKS ON MULTIMAP===============================
//								====================> pathsmmap to arranged_mmap
void arrangement_of_links()
{
	arranged_mmap.clear();
	vector<int>tempvec;

	for (pmmap = pathsmmap.begin(); pmmap != pathsmmap.end(); ++pmmap)
	{
		const vector<int>& vpmmap = pmmap->second;
		//if(!vmatch(ppaths, vapm))				//Condition for eliminate joint links
		{
			for (size_t n = 1; n < vpmmap.size(); ++n)
			{
				for (itgtm = graphtextmap.begin(); itgtm != graphtextmap.end(); ++itgtm)
				{
					const vector<int>& vgtm = itgtm->second;
					{
						if (((vpmmap[n - 1] == vgtm[0]) && (vpmmap[n] == vgtm[1])))
						{
							for (fbwtmp = fiberid_weight_map.begin(); fbwtmp != fiberid_weight_map.end(); ++fbwtmp)
							{
								if (itgtm->first == fbwtmp->first)
								{
									tempvec.push_back(itgtm->first);
								}
							}
						}
					}
				}
			}
		}
		//if (arranged_mmap.size()<MAX_POSSIBLE_PATHS)
		//{
		arranged_mmap.insert(make_pair(pmmap->first, tempvec));
		tempvec.clear();
		//}
		//else
			//break;
	}
}

//DISPLAY WEIGHT_LINKS MULTIMAP
void display_arranged_links()
{
	cout << "\n\nWEIGHT - ARRANGED_LINKS\n";
	cout << "-----------------------";
	for (armmap = arranged_mmap.begin(); armmap != arranged_mmap.end(); ++armmap)
	{
		cout << endl << armmap->first << " => ";
		for (size_t n = 0; n < (*armmap).second.size(); n++)
		{
			cout << " " << armmap->second[n] << " ";
		}
	}
}
//=====================ARRANGING LINKS ON MULTIMAP===============================

//==================START OF FIND BACKUP PATHS FUNCTIONS===============================
bool imatch(const IntVec& vec, int x, int y)
{
	for (size_t i = 1; i < vec.size(); i++)
	{
		if (((vec[i - 1] == x) && (vec[i] == y)) || ((vec[i - 1] == y) && (vec[i] == x)))
		{
			return true;
		}
	}
	return false;
}

//----------------------Find consequence values for map--------------------------------------------
bool vmatch(const IntVec& vec1, const IntVec& vec2)
{
	for (size_t i = 1; i < vec2.size(); i++)
	{
		if (imatch(vec1, vec2[i - 1], vec2[i]))
		{
			return true;
		}
	}
	return false;
}
//===============================END OF FIND BACKUP PATHS FUNCTIONS===========================================

//==============================START CONVERTING NODES TO LINKS========================================
vector<int> convert_to_links(const IntVec& nodevec)
{
	IntVec linkvec;
	linkvec.clear();

	for (size_t n = 1; n < nodevec.size(); ++n)
	{
		for (itgtm = graphtextmap.begin(); itgtm != graphtextmap.end(); ++itgtm)
		{
			const vector<int>& vgtm = itgtm->second;
			{
				if (((nodevec[n - 1] == vgtm[0]) && (nodevec[n] == vgtm[1])))
				{
					//for (fbwtmp = fiberid_weight_map.begin(); fbwtmp != fiberid_weight_map.end(); ++fbwtmp )
					{
						//if(itgtm->first==fbwtmp->first)
						{
							linkvec.push_back(itgtm->first);
						}
					}
				}
			}
		}

	}
	return linkvec;
}
//==============================END CONVERTING NODES TO LINKS========================================

//*****************************PRIMARY SELECTION STARTS***************************************
//------------------------Selecting Primary------------------------------------------------
//SELECTING PRIMARY PATH AND DELETE ALL JOINT PATHS FOR BACKUP
void NEW_select_primary_path()
{
	int flagin = 0;

	cout << "\n\nPRIMARY";
	cout << "\n=========";
	cout << "\nStart...End";
	cout << "\n------------\n";

	for (pmmap = pathsmmap.begin(); pmmap != pathsmmap.end(); ++pmmap)
	{
		ppnodes = pmmap->second;
		pplinks = convert_to_links(ppnodes);	//Converting nodes (ppnodes) to links (pplinks)

		for (int startIndex = 0; startIndex < MAX_SLOTS + 1 - No_Of_Slots; startIndex++)
		{
			if (PrimaryAllAllowRange(pplinks, startIndex, No_Of_Slots, mainmap))
			{
				cout << startIndex << "-" << startIndex + No_Of_Slots - 1 << "  " << "\n";

				prislotposition[Lightpath].push_back(startIndex);						//Inserting lightpath, start position into prislotposition map
				prislotposition[Lightpath].push_back(startIndex + No_Of_Slots - 1);	//Inserting lightpath, end position into prislotposition map
				flagin = 1;
				break;
			}
		}
		if (flagin == 1)
		{
			break;
		}
	}
	if (lpchecker(prislotposition))
	{
		for (pmmap = pathsmmap.begin(); pmmap != pathsmmap.end(); )
		{
			const IntVec& vec = pmmap->second;
			if (vmatch(ppnodes, vec))				//Condition for eliminate joint links
			{
				pathsmmap.erase(pmmap++);
			}
			else
			{
				++pmmap;
			}
		}
	}
	else
	{
		ppnodes.clear();
	}
}

//-----------------------Allowing A Range of Values of ArrayMap if in VecMap--------------------------------
bool PrimaryAllAllowRange(const vector<int>& pplinks, int startIndex, int len, const VecMap& theMap)
{
	for (VecMap::const_iterator it = theMap.begin(); it != theMap.end(); ++it)
	{
		int key = it->first;
		const vector<int>& vec = it->second;

		vector<int>::const_iterator itv = find(pplinks.begin(), pplinks.end(), key);

		if (itv != pplinks.end())
		{
			for (int i = startIndex; i <= startIndex + len - 1; i++)
			{
				if (!Allow(vec[i]))
				{
					return false;
				}
			}
		}
	}
	return true;
}

//----------------------Allowing Special Values and special keys if in VecMap---------------------------------
bool Allow(int value)
{
	if (value == 0)
	{
		return true;
	}
	return false;
}

//********************************PRIMARY SELECTION ENDS***********************************************************


//===========================INSERT INTO SORTED_MAP STARTS================================================================

// returns true if finds value in 'vec'
bool Find(const IntVec& vec, int value)
{
	return find(vec.begin(), vec.end(), value) != vec.end();
}

//---------------TEST A VECTOR NOT IN ANOTHER VECTOR AND IN FIBER_LINK MAP----------------					
bool HaveNoCommonValueEither(const map<int, int>& themap, const IntVec& vec1, const IntVec& vec2)
{
	for (size_t i = 0; i < vec1.size(); i++)
	{
		map<int, int>::const_iterator itthemap = themap.find(vec1[i]);
		int same = itthemap->second;
		{
			if ((Find(vec2, vec1[i])) || (Find(vec2, same)))
			{
				return false;
			}
		}
	}
	return true;
}

//-------------------------------------------------------------------------------------------------
bool AllowKey(const VecMap& primap, int specialKey, int key)
{
	VecMap::const_iterator itpri = primap.find(specialKey);
	const IntVec& specialVec = itpri->second;		//Dereferencing assertion failed error

	itpri = primap.find(key);
	const IntVec& vec = itpri->second;

	if (HaveNoCommonValueEither(fiber_link, specialVec, vec))
	{
		return true;
	}
	return false;
}

//-------------------------------------------------------------------------------------------------------
VecMap CreateNewMap(const VecMap& primarypath, int specialKey, const VecMap& backuppath, const IntVec& specialVec)
{

	sorted_map.clear();
	for (VecMap::const_iterator itsec = backuppath.begin(); itsec != backuppath.end(); ++itsec)
	{
		int key = itsec->first;
		const IntVec& vec = itsec->second;
		if ((key != specialKey) && AllowKey(primarypath, specialKey, key))
		{
			for (size_t n = 0; n < vec.size(); n++)
			{
				if (Find(specialVec, (vec[n])))
				{
					sorted_map[vec[n]].push_back(key);
				}
			}
		}
	}
	return sorted_map;
}

//-----------------------------------DISPLAYING SORTED_MAP------------------------------------
void display_sorted_map()
{
	cout << "\n\nSORTED_MAP CONTAINS POSSIBLE SHARABLE LIGHTPATHS(LINKS=>LIGHTPATH)\n";
	cout << "-------------------------------------------------------------------";
	for (VecMap::const_iterator amap = sorted_map.begin(); amap != sorted_map.end(); ++amap)
	{
		cout << endl << amap->first << "=>";
		for (size_t n = 0; n < (*amap).second.size(); n++)
		{
			cout << " " << amap->second[n] << " ";
		}
	}
	cout << "\n\n";
}
//----------------------------------------------------------------------------------------------
//==============================INSERT INTO SORTED_MAP ENDS========================================================

//======================Inserting other shaerable id inside splvec========================================

//-------------------Other Shareable ids from shared_id map and whether in sorted_map values-------------
bool IsValueInMap(int value, const VecMap& m2)
{
	for (VecMap::const_iterator itm = m2.begin(); itm != m2.end(); ++itm)
	{
		const IntVec& vec = itm->second;
		IntVec::const_iterator itmv = find(vec.begin(), vec.end(), value);
		if (itmv != vec.end())
		{
			return true;
		}
	}
	return false;
}

bool AreAllValuesInMap(const IntVec& v, const VecMap& m2)
{
	for (size_t i = 0; i < v.size(); i++)
	{
		if (!IsValueInMap(v[i], m2))
		{
			return false;
		}
	}
	return true;
}
//-------------------Other Shareable ids from shared_id map and whether in sorted_map values-------------

IntVec other_shareable_ids(const VecMap& sorted_map)
{
	splvec.clear();
	for (slpsid = shared_id.begin(); slpsid != shared_id.end(); ++slpsid)
	{
		const IntVec& myv1 = slpsid->second;
		if (AreAllValuesInMap(myv1, sorted_map))
		{
			splvec.push_back(slpsid->first);
		}
	}
	return splvec;
}
//======================Inserting other shaerable id inside splvec======================================

//==================================BACKUP SELECTION WITH SHARING STARTS======================================
//-----------------------Allowing A Range of Values of ArrayMap if in VecMap--------------------------------
bool SecondaryAllAllowRange(const IntVec& mvec, int startIndex, int len, const VecMap& theMap, const VecMap& m1, const IntVec& v)
{

	for (VecMap::const_iterator it = theMap.begin(); it != theMap.end(); ++it)
	{
		int key = it->first;
		const IntVec& vec = it->second;

		IntVec::const_iterator itv = find(mvec.begin(), mvec.end(), key);
		if (itv != mvec.end())
		{

			for (int i = startIndex; i <= startIndex + len - 1; i++)
			{
				if (!SecondaryAllow(vec[i], key, m1, v))
				{
					return false;
				}
			}
		}
	}
	return true;
}

//----------------------Allowing Special Values if in VecMap---------------------------------
bool SecondaryAllow(int value, int key, const VecMap& m1, const IntVec& v)
{
	if (value == 0)
	{
		return true;
	}
	IntVec::const_iterator itv = find(v.begin(), v.end(), value);
	if (itv != v.end())
	{
		return true;
	}
	VecMap::const_iterator itm = m1.find(key);
	if (itm != m1.end())
	{
		const IntVec& vec = itm->second;
		return (find(vec.begin(), vec.end(), value) != vec.end());
	}
	return false;
}

//-----------------Counting Special Values of ArrayMap for Certain Range if in VecMap------------------
int SecondaryCountAllSpecialInRange(const IntVec& mvec, int startIndex, int len, VecMap& theMap, const VecMap& m1, const IntVec& v)
{
	int count = 0;
	for (VecMap::const_iterator it = theMap.begin(); it != theMap.end(); ++it)
	{
		int key = it->first;
		const IntVec& vec = it->second;

		IntVec::const_iterator itv = find(mvec.begin(), mvec.end(), key);
		if (itv != mvec.end())
		{

			for (int i = startIndex; i <= startIndex + len - 1; i++)
			{
				if (SecondarySpecial(vec[i], m1, v))
					count++;
			}
		}
	}
	return count;
}

//------------------------Allowing Special Values if in VecMap-----------------------------------
bool SecondarySpecial(int value, const VecMap& m1, const IntVec& v)
{
	IntVec::const_iterator itv = find(v.begin(), v.end(), value);
	if (itv != v.end())
	{
		return true;
	}
	for (VecMap::const_iterator itm = m1.begin(); itm != m1.end(); ++itm)
	{
		const IntVec& myvec = itm->second;
		IntVec::const_iterator itv = find(myvec.begin(), myvec.end(), value);
		if (itv != myvec.end())
		{
			return true;
		}
	}
	return false;
}

//--------------------------------------------------------------------------------------------------
VecMul CreateNewVector(const IntVec& myvec, VecMap& mainmap, const VecMap& newmap, const IntVec& vect)
{
	VMap.clear();
	IntVec IVec;

	int countMax = 0;
	int indexMax = -1;
	int Entry = 1;		//For counting the possible entries

	cout << "\nSECONDARY";
	cout << "\n=========";
	cout << "\nEntry -> Start-End -> Count";
	cout << "\n---------------------------";

	for (int startIndex = 0; startIndex < MAX_SLOTS + 1 - No_Of_Slots; startIndex++)  // 11 + 1 - 3 = 9 => 8-10
	{
		if (SecondaryAllAllowRange(myvec, startIndex, No_Of_Slots, mainmap, newmap, vect))
		{
			int count = SecondaryCountAllSpecialInRange(myvec, startIndex, No_Of_Slots, mainmap, newmap, vect);

			cout << "\n" << Entry++ << "->" << startIndex << "-" << startIndex + No_Of_Slots - 1 << "->" << count;

			IVec.clear();
			IVec.push_back(startIndex);
			IVec.push_back(startIndex + No_Of_Slots - 1);

			VMap.insert(make_pair(count, IVec));

		}
	}
	return VMap;
}

//==================================BACKUP SELECTION WITH SHARING ENDS========================================

//**************************INSERTING LIGHTPATHS INTO MAINMAP STARTS******************************************
//-----------------------------------------------------------------------------
int FindMapKeyWithSameValues(VecMap& smap, IntVec& vec)
{
	for (VecMap::iterator its = smap.begin(); its != smap.end(); ++its)
	{
		int key = its->first;
		IntVec& sVec = its->second;
		sort(sVec.begin(), sVec.end());
		sort(vec.begin(), vec.end());
		if (sVec == vec)
		{
			return key;
		}
	}
	return -1;
}

//--------------------------------------------------------------------------------
void ProcessData(IntVec& mVec, VecMap& smap, int j, int positionKey)
{
	int val = mVec[j];

	if (val == 0)
	{
		mVec[j] = positionKey;
		return;
	}

	VecMap::iterator its = smap.find(val); // find key equal to val
	if (its != smap.end())
	{
		int skey = its->first;
		IntVec& sVec = its->second; // note reference
		IntVec newVec;
		newVec = sVec;

		if (!Find(newVec, positionKey))
		{
			newVec.push_back(positionKey);
			int key = FindMapKeyWithSameValues(smap, newVec);
			if (key < 0) // not found
			{
				VecMap::iterator its = smap.end();
				its--;
				key = its->first + 1; // new key
				smap[key] = newVec;
			}
			mVec[j] = key;
			return;
		}

		mVec[j] = skey;
		return;
	}

	if (val < SHLPIDSTART)
	{
		IntVec vec; // new vector with two entries
		vec.push_back(val);
		vec.push_back(positionKey);
		int key = FindMapKeyWithSameValues(smap, vec);
		if (key < 0) // not found
		{
			if (smap.empty())
			{
				key = SHLPIDSTART;
			}
			else
			{
				VecMap::iterator its = smap.end();
				its--;
				key = its->first + 1; // new key
			}
			smap[key] = vec;
		}
		mVec[j] = key;
		return;
	}
}

//---------------------------------Filling backup paths in mainmap----------------------------------------
void insert_backup_lp(const IntVec& backvec)
{
	bsltposn = backslotposition.find(Lightpath);
	if (bsltposn != backslotposition.end())
	{
		IntVec& posvec = bsltposn->second;
		int start = posvec[0];
		int end = posvec[1];
		for (itmm = mainmap.begin(); itmm != mainmap.end(); ++itmm)
		{
			int itmmkey = itmm->first;
			IntVec& itmmvec = itmm->second;
			IntVec::const_iterator itv = find(backvec.begin(), backvec.end(), itmmkey);
			if (itv != backvec.end())
			{
				for (int index = start; index <= end; index++)
				{
					ProcessData(itmmvec, shared_id, index, Lightpath);			//int itmmkey is only for splitMap operation	
				}
			}
		}
	}
}

//---------------------------------Filling primary paths in mainmap----------------------------------------
void insert_primary_lp(IntVec& primvec)
{
	psltposn = prislotposition.find(Lightpath);
	if (psltposn != prislotposition.end())
	{
		IntVec& posvec = psltposn->second;
		int start = posvec[0];
		int end = posvec[1];
		for (itmm = mainmap.begin(); itmm != mainmap.end(); ++itmm)
		{
			int itmmkey = itmm->first;
			IntVec& itmmvec = itmm->second;
			IntVec::const_iterator itv = find(primvec.begin(), primvec.end(), itmmkey);
			if (itv != primvec.end())
			{
				fill(itmmvec.begin() + start, itmmvec.begin() + start + No_Of_Slots, Lightpath);
			}
		}
	}
	//pplinks.clear();
}


//********************************SECONDARY SELECTION STARTS***************************************
void select_slots_for_secondarypath()
{
	VecMap sorted_map;
	IntVec vec1;

	if (!arranged_mmap.empty())
	{
		primarypath[Lightpath] = pplinks;			//Insert lightpaths=>links into primarypath map

		const IntVec& specialVec = arranged_mmap.begin()->second;

		sorted_map = CreateNewMap(primarypath, Lightpath, backuppath, specialVec);   //specialVec is a selected backup paths from wtmap

		display_sorted_map();			//Displaying sorted_map

		vec1 = other_shareable_ids(sorted_map);			//Inserting other shareable ids from shared_id map into splvec

		cout << "\nSPECIAL VALUES TO BE SHARED";
		cout << "\n---------------------------\n";
		for (vector<int>::iterator itve = splvec.begin(); itve != splvec.end(); ++itve)
		{
			cout << (*itve) << " ";
		}
		cout << "\n\n";

		VMap = CreateNewVector(specialVec, mainmap, sorted_map, vec1);			//specialVec is a selected backup paths

		display_slots_range();		//Displaying count and slots range.

		/*
			cout<<"\nSELECTED SECONDARY";
			cout<<"\n==================";
			cout<<"\nStart...End...Max";
			cout<<"\n-----------------\n";

			for(IntVec::const_iterator itv=vec.begin(); itv!=vec.end(); ++itv)
			{
				cout<<(*itv)<<" ";
			}
			cout<<"\n\n";
		*/

		/*
		if(!VMap.empty())
		{
			IntVec slotVec;
			if(VMap.rbegin()->first!=0)
			{
				slotVec=VMap.rbegin()->second;
			}
			else if(VMap.rbegin()->first==0)
			{
				slotVec=VMap.begin()->second;
			}

			cout<<"\n";
			for(IntVec::iterator revIT=slotVec.begin(); revIT!=slotVec.end();++revIT)
			{
				cout<<(*revIT)<<" ";
			}

			NodeLayer=ProcessLayer(slotVec[0], slotVec[1]);			//For defining layers in nodes (implementing Flex-Wss)

			if(NodeLayer.size()==bpnodes.size())			//For defining layers in nodes (implementing Flex-Wss)
			{
				NotifySet1=PBPS_COND_1(slotVec[0], slotVec[1]);			//For the ADDITIONAL CONDITION (implementing PBPS)
				NotifySet2=PBPS_COND_2(slotVec[0], slotVec[1]);			//For the ADDITIONAL CONDITION (implementing PBPS)

				if((NotifySet1.size()==bplinks.size()-1)&&(NotifySet2.size()==bplinks.size()-1))
				{
					insert_lightpath_into_layer();			//Insert nodes and layers into mainnodemap and layer_id

					backslotposition[Lightpath].push_back(slotVec[0]);		//Inserting lightpath, start position into backslotposition map
					backslotposition[Lightpath].push_back(slotVec[1]);		//Inserting lightpath, end position into backslotposition map
				}
			}
		}
		*/

		//================FOR CHECKING ALL AVAILABLE RANGES=============================================

		if (!VMap.empty())
		{
			int timer = 0;
			int timerACon = 0;
			if (VMap.rbegin()->first != 0)
			{
				for (VecMul::reverse_iterator revIT = VMap.rbegin(); revIT != VMap.rend(); ++revIT)
				{
					IntVec& slotVec = revIT->second;

					NodeLayer = ProcessLayer(slotVec[0], slotVec[1]);			//For defining layers in nodes (implementing Flex-Wss)

					if (NodeLayer.size() == bpnodes.size())			//For defining layers in nodes (implementing Flex-Wss)
					{
						NotifySet1 = PBPS_COND_1(slotVec[0], slotVec[1]);			//For the ADDITIONAL CONDITION (implementing PBPS)
						NotifySet2 = PBPS_COND_2(slotVec[0], slotVec[1]);			//For the ADDITIONAL CONDITION (implementing PBPS)

						if ((NotifySet1.size() == bplinks.size() - 1) && (NotifySet2.size() == bplinks.size() - 1))
						{
							//for(set<int>::iterator iT = NotifySet.begin(); iT != NotifySet.end(); ++iT)
							//{
							//if(*iT==1)

							insert_lightpath_into_layer();			//Insert nodes and layers into mainnodemap and layer_id

							backslotposition[Lightpath].push_back(slotVec[0]);		//Inserting lightpath, start position into backslotposition map
							backslotposition[Lightpath].push_back(slotVec[1]);		//Inserting lightpath, end position into backslotposition map
							break;
						}
						else
						{
							timerACon++;
						}
					}
					else
					{
						timer++;
					}
				}
			}
			else if (VMap.rbegin()->first == 0)
			{
				for (VecMul::iterator fwdIT = VMap.begin(); fwdIT != VMap.end(); ++fwdIT)
				{
					IntVec& slotVec = fwdIT->second;

					NodeLayer = ProcessLayer(slotVec[0], slotVec[1]);			//For defining layers in nodes (implementing Flex-Wss)

					if (NodeLayer.size() == bpnodes.size())			//For defining layers in nodes (implementing Flex-Wss)
					{
						NotifySet1 = PBPS_COND_1(slotVec[0], slotVec[1]);			//For the ADDITIONAL CONDITION (implementing PBPS)
						NotifySet2 = PBPS_COND_2(slotVec[0], slotVec[1]);			//For the ADDITIONAL CONDITION (implementing PBPS)

						if ((NotifySet1.size() == bplinks.size() - 1) && (NotifySet2.size() == bplinks.size() - 1))
						{
							//for(set<int>::iterator iT = NotifySet.begin(); iT != NotifySet.end(); ++iT)
							//{
							//if(*iT==1)

							insert_lightpath_into_layer();			//Insert nodes and layers into mainnodemap and layer_id

							backslotposition[Lightpath].push_back(slotVec[0]);		//Inserting lightpath, start position into backslotposition map
							backslotposition[Lightpath].push_back(slotVec[1]);		//Inserting lightpath, end position into backslotposition map
							break;
						}
						else
						{
							timerACon++;
						}
					}
					else
					{
						timer++;
					}
				}
			}
			if (timer == VMap.size())
			{
				NOLAYERV.push_back(Lightpath);
			}
			if (timerACon == VMap.size())
			{
				NOADDCONV.push_back(Lightpath);
			}
			if ((timer != 0) && (timerACon != 0) && (timer + timerACon == VMap.size()))
			{
				XBLOCKV.push_back(Lightpath);	//Block either no layer/no add. con. satisfied
			}
		}
		else if (VMap.empty())
		{
			NOBACKSLOTV.push_back(Lightpath);
		}
		//================FOR CHECKING ALL AVAILABLE RANGES=============================================

	}
	display_node_layer();
	cout << "\n";
}

//----------------------------------------------------------------------------------------
void display_slots_range()
{
	cout << "\n\n\nSLOTS RANGE CONTAINS (COUNT=>SLOTS)\n";
	cout << "-----------------------------------------";
	for (VecMul::iterator it = VMap.begin(); it != VMap.end(); ++it)		//Display shared_id and lightpaths
	{
		cout << endl << (*it).first << "=>";
		for (size_t n = 0; n < (*it).second.size(); n++)
		{
			cout << " " << it->second[n] << " ";
		}
	}
}
//********************************SECONDARY SELECTION ENDS***************************************

//***********************WHOLE INSERTION OF PRIMARY AND BACKUP PATHS STARTS******************************
void insert_all_lightpaths()
{
	insert_primary_lp(pplinks);				//Insert primary slots in mainmap

	insert_backup_lp(bplinks);				//Insert backup slots in mainmap

	//primarypath[Lightpath] = pplinks;			//Insert lightpaths=>links into primarypath map		//ALREADY INSERTED

	primarynode[Lightpath] = ppnodes;			//Insert lightpaths=>nodes into primarynode map

	backuppath[Lightpath] = bplinks;		//Insert lightpaths=>links into backuppath map

	backupnode[Lightpath] = bpnodes;		//Insert lightpaths=>nodes into backupnode map

	//-------------------------FOR LAYERS (FLEX-WSS)-----------------------------------------
	insert_lightpath_into_node_lps();		//Insert nodes and lightpaths in node_lps map
	//-------------------------FOR LAYERS (FLEX-WSS)-----------------------------------------

	//----------------------------FOR PBPS------------------------------------------------------
	insert_splitmulmap(bplinks);			//Insert splitmap for link=>cur_lp->split_lps details

	insert_lightpath_into_link_lps();		//Inserting bplinks and lightpaths in bplink_lps map
	//----------------------------FOR PBPS-------------------------------------------------------

}
//**********************WHOLE INSERTION OF PRIMARY AND BACKUP PATHS ENDS*********************************

//**************************DELETION FUNCTION STARTS**********************************************************
//===========================FOR PRIMARY DELETION=============================================================
//DELETE LIGHTPATH IN LINKS FROM MAINMAP MULTIMAP
void del_pri_lp_from_mainmap(int dLightpath)
{
	int positionKey = dLightpath;

	VecMap::iterator pit = primarypath.find(dLightpath);
	if (pit != primarypath.end())
	{
		IntVec& vect = pit->second;
		VecMap::iterator pos = prislotposition.find(positionKey);
		if (pos != prislotposition.end())
		{
			IntVec& posvec = pos->second;
			int start = posvec[0];
			int end = posvec[1];
			for (VecMap::iterator itmm = mainmap.begin(); itmm != mainmap.end(); ++itmm)
			{
				int key = itmm->first;
				IntVec& vec = itmm->second;
				IntVec::iterator itv = find(vect.begin(), vect.end(), key);
				if (itv != vect.end())
				{
					for (int j = start; j <= end; j++)
					{
						ProcessDeleteData(vec, j, positionKey);
					}
				}
			}
		}
	}
}

void ProcessDeleteData(IntVec& mVec, int j, int positionKey)
{
	int val = mVec[j];

	if (val == positionKey)
	{
		mVec[j] = 0;
		return;
	}
}

//===========================FOR PRIMARY DELETION===================================

//(((((((((((((((-------------START FOR BACKUP WITH SHARING DELETION------------))))))))))))))))))))))))))))
//===========================DELETE SHARING LPS MAIN FUNCTION STARTS========================================
void del_back_lp_from_mainmap(int dLightpath)
{

	int value = dLightpath;

	map<int, int> imap; // map to hold substitution values

	imap[value] = 0;

	VecMap::iterator iback = backuppath.find(dLightpath);
	if (iback != backuppath.end())
	{
		IntVec& dellinks = iback->second;


		if (!dellinks.empty())		//Checks if dLightpath available in mypathbackup
		{
			del_lightpath_from_bplink_lps(dellinks, dLightpath);	//For deleting lightpaths from bplink_lps map (FOR PBPS)

			for (map<int, vector<int>>::const_iterator ip = backslotposition.begin(); ip != backslotposition.end();) // map-erase idiom coming
			{
				int value = ip->first;
				const vector<int>& vec = ip->second;
				assert(vec.size() == 2);
				pair<int, int> position = make_pair(vec[0], vec[1]);
				if (value == dLightpath)
				{
					backslotposition.erase(ip++);

					for (map<int, vector<int>>::iterator slpsid = shared_id.begin(); slpsid != shared_id.end(); ) // map-erase idiom coming
					{
						int key = slpsid->first;
						vector<int>& vec = slpsid->second;
						sort(vec.begin(), vec.end()); // sort so can check same elements using operator ==()

						if (Find(vec, value)) // 'value' is in 'vec'
						{
							vec.erase(remove(vec.begin(), vec.end(), value), vec.end());  // remove 'value' from vector
							if (vec.size() == 1) // modified vector has only one entry
							{
								if (Find(dellinks, mainmap, key, position))	//'key' is in 'mmap'
								{
									imap[key] = vec[0]; // substitution value
								}
								shared_id.erase(slpsid++);   // remove vector from 'smap'
							}
							else
							{
								if (Find(dellinks, mainmap, key, position))
								{
									imap[key] = NewValue(key, vec, shared_id); // substitution value
								}
								++slpsid;
							}
						}
						else
						{
							++slpsid;
						}
					}
				} //end of if (value == dLightpath)
				else
				{
					++ip;
				}

				for (map<int, int>::iterator it = imap.begin(); it != imap.end(); ++it)
				{
					Replace(dellinks, mainmap, it->first, it->second, position);
				}
			}
			dellinks.clear();

			DelVec.push_back(dLightpath);	//Insert deleting lightpath into a DelVec vector
		} //end of if(!dellinks.empty())
	} //end of if (iback != backuppath.end())
	else
	{
		cout << "\n\nNO LIGHTPATH EXIST TO BE DELETED!!!!\n";
	}

	//--------------------Displaying replacing values after deleting a lightpath--------------------------
	cout << "\n\nSUBSTITUTE MAP (FOR STORING REPLACEMENT VALUES AFTER DELETING A LIGHTPATHS)\n";
	cout << "-----------------------------------------------------------------------------\n";
	for (map<int, int>::iterator it = imap.begin(); it != imap.end(); ++it)
	{
		cout << it->first << " => " << it->second << endl;
	}
	cout << "\n\n";

}
//------------------------------------------------------------------------------------------------------
//===========================DELETE SHARING LPS MAIN FUNCTION ENDS========================================

//--------------------------FOR DELETE SHARING LPS SUB FUNCTION-----------------------------------------------
//-----------------------------------------------------------------------------
bool Find(const vector<int>& myvector, const map<int, vector<int>>& vecMap, int value, const pair<int, int>& position)
{
	for (map<int, vector<int>>::const_iterator it = vecMap.begin(); it != vecMap.end(); ++it)
	{
		int key = it->first;
		const vector<int>& vec = it->second;
		vector<int>::const_iterator itv = find(myvector.begin(), myvector.end(), key);
		if (itv != myvector.end())
		{
			vector<int>::const_iterator beg = vec.begin() + position.first;
			vector<int>::const_iterator end = vec.begin() + position.second + 1;
			if (find(beg, end, value) != end)
			{
				return true;
			}
		}
	}
	return false;
}
//----------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
// returns different key of smap if vectors match, otherwise return 'key'
// both vec and vectors in smap are assumed to be sorted
int NewValue(int key, const vector<int>& vec, const map<int, vector<int>>& smap)
{
	for (map<int, vector<int>>::const_iterator it1 = smap.begin(); it1 != smap.end(); ++it1)
	{
		int key1 = it1->first;
		const vector<int>& vec1 = it1->second;
		if ((key1 != key) && (vec1 == vec))
		{
			return key1;
		}
	}
	return key;
}
//-----------------------------------------------------------------------------------
//-------------------------------------------------------------------------------

void Replace(const vector<int>& dellinks, map<int, vector<int>>& mainmap, int oldValue, int newValue, const pair<int, int>& position)
{
	for (map<int, vector<int>>::iterator it = mainmap.begin(); it != mainmap.end(); ++it)
	{
		int key = it->first;
		vector<int>& vec = it->second;
		vector<int>::const_iterator itv = find(dellinks.begin(), dellinks.end(), key);
		if (itv != dellinks.end())
		{
			vector<int>::iterator beg = vec.begin() + position.first;
			vector<int>::iterator end = vec.begin() + position.second + 1;
			//if (Find(mmap, oldValue, position))		
			replace(beg, end, oldValue, newValue);
		}
	}
}
//--------------------------------------------------------------------------------
//----------------------------FOR DELETE SHARING LPS SUB FUNCTION------------------------------------------

//DISPLAY DELETED LIGHTPATHS
void display_deleted_lightpaths()
{
	cout << "\nDELETED LIGHTPATHS\n";
	cout << "--------------------\n";
	for (IntVec::iterator deVec = DelVec.begin(); deVec != DelVec.end(); ++deVec)
	{
		cout << (*deVec) << ", ";
	}
	cout << "\n\n";
}

//((((((((((((((((((((-------------END SHARING DELETION FUNCTION------------))))))))))))))))))))))))))))))))))

//*********************************DELETION FUNCTION ENDS****************************************************

//===========================DISPLAYING FUNCTIONS START===================================
//DISPLAY SHARED_ID MAP
void display_shared_id()
{
	cout << "\n\n\nSHARED ID CONTAINS (SHARED-ID=>LIGHTPATHS)\n";
	cout << "-----------------------------------------";
	for (slpsid = shared_id.begin(); slpsid != shared_id.end(); ++slpsid)		//Display shared_id and lightpaths
	{
		cout << endl << (*slpsid).first << "=>";
		for (size_t n = 0; n < (*slpsid).second.size(); n++)
		{
			cout << " " << slpsid->second[n] << " ";
		}
	}
}

//DISPLAY PSLOTPOSITION MAP
void display_pslotposition()
{
	cout << "\n\n\nPRIMARY SLOTS POSITIONS (PRIMARY LIGHTPATH=>FROM-TO)\n";
	cout << "----------------------------------------------------";
	for (psltposn = prislotposition.begin(); psltposn != prislotposition.end(); ++psltposn)
	{
		cout << endl << psltposn->first << "=>";
		for (size_t n = 0; n < (*psltposn).second.size(); n++)
		{
			cout << " " << psltposn->second[n] << " ";
		}
	}
}

//DISPLAY BSLOTPOSITION MAP
void display_bslotposition()
{
	cout << "\n\n\nBACKUP SLOTS POSITIONS (BACKUP LIGHTPATH=>FROM-TO)\n";
	cout << "--------------------------------------------------";
	for (bsltposn = backslotposition.begin(); bsltposn != backslotposition.end(); ++bsltposn)
	{
		cout << endl << bsltposn->first << "=>";
		for (size_t n = 0; n < (*bsltposn).second.size(); n++)
		{
			cout << " " << bsltposn->second[n] << " ";
		}
	}
}

//DISPLAY PRIMARYPATH MAP
void display_primarypath()
{
	cout << "\n\n\nPRIMARY PATH CONTAINS (PRIMARY LIGHTPATH=>LINKS)\n";
	cout << "-----------------------------------------------";
	for (itppath = primarypath.begin(); itppath != primarypath.end(); ++itppath)		//Display primary paths and links of map primarypath
	{
		cout << endl << (*itppath).first << "=>";
		for (size_t n = 0; n < (*itppath).second.size(); n++)
		{
			cout << " " << itppath->second[n] << " ";
		}
	}
}

//DISPLAY PRIMARYNODE MAP
void display_primarynode()
{
	cout << "\n\n\nPRIMARY NODE CONTAINS (PRIMARY LIGHTPATH=>NODES)\n";
	cout << "-----------------------------------------------";
	for (itpnode = primarynode.begin(); itpnode != primarynode.end(); ++itpnode)		//Display primary paths and nodes of map primarynode
	{
		cout << endl << (*itpnode).first << "=>";
		for (size_t n = 0; n < (*itpnode).second.size(); n++)
		{
			cout << " " << itpnode->second[n] << " ";
		}
	}
}

//DISPLAY BACKUPPATH MAP
void display_backuppath()
{
	cout << "\n\n\nBACKUP PATH CONTAINS (BACKUP LIGHTPATH=>LINKS)\n";
	cout << "----------------------------------------------";
	for (itbpath = backuppath.begin(); itbpath != backuppath.end(); ++itbpath)		//Display backup paths and links of map backuppath
	{
		cout << "\n" << (*itbpath).first << "=>";
		for (size_t n = 0; n < (*itbpath).second.size(); n++)
		{
			cout << " " << itbpath->second[n] << " ";
		}
	}
}

//DISPLAY BACKUPNODE MAP
void display_backupnode()
{
	cout << "\n\n\nBACKUP NODE CONTAINS (BACKUP LIGHTPATH=>NODES)\n";
	cout << "----------------------------------------------";
	for (itbnode = backupnode.begin(); itbnode != backupnode.end(); ++itbnode)		//Display backup paths and nodes of map backupnode
	{
		cout << "\n" << (*itbnode).first << "=>";
		for (size_t n = 0; n < (*itbnode).second.size(); n++)
		{
			cout << " " << itbnode->second[n] << " ";
		}
	}
}

//DISPLAY PRIMARY LINKS												
void display_pplinks()
{
	cout << "\n\nPRIMARY LINKS\n";
	cout << "---------------\n";
	for (IntVec::iterator it = pplinks.begin(); it != pplinks.end(); ++it)
	{
		cout << " " << (*it);
	}
}
//===========================DISPLAYING FUNCTIONS END===================================

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---FOR TRIAL STARTS---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//===========================FOR ADMITTING A LIGHTPATH===================================
void admit_lightpath(int start, int end, int Lightpath)

//if((No_Of_Slots!=0)&&(start!=0)&&(end!=0)&&(start!=end))
{

	findallpaths(start, end);

	display_all_paths_mmap();

	collection_of_paths();

	display_pathmaps();

	//Before deleting joint paths
	if (pathsmmap.empty())
	{
		RATESBLOCKED.push_back(Data_Rate);
		NOPRIPATH.insert(Lightpath);
		printcase1(); //NO PRIMARY PATH AVAILABLE (CASE-1)
	}

	//============================================================================
	NEW_select_primary_path();

	display_pathmaps();			//Displaying only backup paths
//============================================================================

	arrangement_of_links();		//storing links according to the weight from  pathsmap map

	display_arranged_links();


	if (ppnodes.empty())
	{
		RATESBLOCKED.push_back(Data_Rate);
		NOPRISLOT.insert(Lightpath);
		//printcase2();	//NO SLOTS FOR PRIMARY PATHS
	}

	else if (pathsmmap.empty())
	{
		RATESBLOCKED.push_back(Data_Rate);
		NOBACKPATH.insert(Lightpath);
		printcase3(); //NO BACKUP PATH AVAILABLE (CASE-3)
	}

	else if (pathsmmap.size() >= 1)
	{

		int countp = 0;

		//cout<<"\n\nSHUUUUUUUUUUUUUUUUUUTHA\n  ";
		//for(IntVec::iterator ite = ppnodes.begin(); ite != ppnodes.end(); ++ite)
		//{
			//cout<<(*ite)<<" ";
		//}
		//cout<<" SHEEEEEEEEEEEEEEEEEEEL\n";

		if (lpchecker(prislotposition))
		{
			display_pplinks();

			int noslot = 0;
			int nopath = 0;

			NOBACKSLOTV.clear();
			NOLAYERV.clear();
			NOADDCONV.clear();
			XBLOCKV.clear();

			for (VecMul::iterator pmmap1 = pathsmmap.begin(); pmmap1 != pathsmmap.end(); ++pmmap1)
			{
				bpnodes = pmmap1->second;

				cout << "\n\nBACKUP - " << pmmap1->first;
				bplinks = convert_to_links(bpnodes);	//Converting nodes (bpnodes) to links (bplinks)
				select_slots_for_secondarypath();

				if (lpchecker(backslotposition))
				{
					cout << "BACKUP - " << pmmap1->first << "  ";
					insert_all_lightpaths();
					//flag=1; 
					break;
				}
				else // if(!lpchecker(backslotposition))
				{
					nopath++;
				}
			}

			/*
			if (nopath==(pathsmmap.size()))
			{
				RATESBLOCKED.push_back(Data_Rate);
				printcase3();	//NO BACKUP PATHS AVAILABLE
				NOBACKPATH.insert(Lightpath);
				prislotposition.erase(Lightpath);
				//flag=1;
				//break;
			}
			*/

			//if((noslot==pathsmmap.size()-1))
			//{			
				//primarypath.erase(Lightpath);
				//printcase4();	//NO SLOTS FOR BACKUP PATHS
				//break;
			//}

			//if((noslot!=0)&&(nopath!=0)&&((noslot+nopath)==pathsmmap.size()))
			//{
				//RATESBLOCKED.push_back(Data_Rate);
				//XBLOCK.insert(Lightpath);
				//prislotposition.erase(Lightpath);
				//cout<<"\n\nEither no backup/no slots for backup";
				//flag=1;
				//break;
			//}

			//-----(1) For only "no slots for backup path"
			if (NOBACKSLOTV.size() == pathsmmap.size())
			{
				RATESBLOCKED.push_back(Data_Rate);
				printcase4();	//NO SLOTS FOR BACKUP PATHS
				NOBACKSLOT.insert(Lightpath);
				prislotposition.erase(Lightpath);
				//flag=1;
				//break;
			}

			//-----(2) For only "no layers"
			if (NOLAYERV.size() == pathsmmap.size())
			{
				RATESBLOCKED.push_back(Data_Rate);
				printcase5();	//NO LAYERS AVAILABLE
				NOLAYER.insert(Lightpath);
				prislotposition.erase(Lightpath);
				//flag=1;
				//break;
			}

			//-----(3) For only "no additional condition satisfied"
			if (NOADDCONV.size() == pathsmmap.size())
			{
				RATESBLOCKED.push_back(Data_Rate);
				printcase6();	//NO ADDITIONAL CONDITION SATISFIED
				NOADDCON.insert(Lightpath);
				prislotposition.erase(Lightpath);
				//flag=1;
				//break;
			}

			//-----For either three of the above
			if (XBLOCKV.size() == pathsmmap.size())
			{
				RATESBLOCKED.push_back(Data_Rate);
				//printcase5();	//NO LAYERS AVAILABLE
				XBLOCK.insert(Lightpath);
				prislotposition.erase(Lightpath);
				//flag=1;
				//break;
			}
		}
	}

	cout << "\n\n";

	cout << "\n\n*******************************************************************";
	cout << "\n----------------------AFTER ADDING A LIGHTPATH---------------------";
	cout << "\n*******************************************************************";
	display_primarynode();		//display primarynode map
	display_primarypath();		//display primarypath map
	display_pslotposition();	//display slot positions of primary lightpath
	display_backupnode();		//display backupnode map
	display_backuppath();		//display backuppath map
	display_bslotposition();	//display slot positions of backup lightpath
	display_shared_id();	//display sharing details of shared_id
	display_mainmap();		//display the whole lightpath and their slots which are in mainmap map

	displaying_layer_id();		//displaying the whole id for mainnodemap and their lightpaths for each layer
	//display_mainnodemap();		//display the whole nodes and their layers which are in mainnodemap map
	display_node_lps();		//Displaying nodes and their current lightpath

	displaying_bplink_lps();	//displaying bplink and lightpaths
	display_splitmulmap();	//display link->lightpaths->splitting_lps details

}
//===========================FOR ADMITTING A LIGHTPATH===================================

//===========================FOR REMOVING A LIGHTPATH===================================
void depart_lightpath(int dLightpath)
{

	//--------------Automated Deletion----------------------------------
		//if(backuppath.size()>MAX_TRAFFIC)
		//{
			//VecMap::iterator itbsp = backslotposition.begin();
			//advance(itbsp, rand() % backslotposition.size());
			//int random_key = itbsp->first;
			//dLightpath=random_key;

			//dLightpath = (rand() % (backslotposition.end()->first- backslotposition.begin->first + 1 ) + backslotposition.begin->first);
			//dLightpath=rand() % 100;
	cout << "\n\nLightpath to be removed: " << dLightpath;
	//--------------Automated Deletion----------------------------------

			//cout<<"\nEnter the Delete Lightpath: ";
			//cin>>dLightpath;

	//------------------------FOR PRIMARY DELETION STARTS-----------------------------------
	del_pri_lp_from_mainmap(dLightpath);	//delete primary lightpath from mainmap map 
	primarynode.erase(dLightpath);		//delete primary nodes from primarynode map
	primarypath.erase(dLightpath);;		//delete primary lightpath from primarypath map
	prislotposition.erase(dLightpath);	//delete primary lightpath from prislotposition map
	//------------------------FOR PRIMARY DELETION ENDS-----------------------------------

	//------------------------FOR DELETING LIGHTPATH IN LAYERS STARTS-----------------------------------
	del_lightpath_in_layers(dLightpath);
	//------------------------FOR DELETING LIGHTPATH IN LAYERS ENDS-----------------------------------
	//------------------------FOR NODE_LPS DELETION STARTS-----------------------------------
	del_lightpath_in_node_lps(dLightpath);
	//------------------------FOR NODE_LPS DELETION ENDS-----------------------------------

	//------------------------FOR SHARED BACKUP DELETION STARTS-----------------------------------
	del_back_lp_from_mainmap(dLightpath);	//delete backup lightpath from mainmap map, (backslotposition map delete included)
	backupnode.erase(dLightpath);		//delete backup nodes from backupnode map
	backuppath.erase(dLightpath);;		//delete backup lightpath from backuppath map
	//------------------------FOR SHARED BACKUP DELETION ENDS-----------------------------------

	//------------------------FOR SPLITMULMAP DELETION STARTS-----------------------------------
	del_lightpath_from_splitmulmap(dLightpath);
	//------------------------FOR SPLITMULMAP DELETION STARTS-----------------------------------

	cout << "\n\n*********************************************************************";
	cout << "\n----------------------AFTER DELETING A LIGHTPATH---------------------";
	cout << "\n*********************************************************************";
	display_primarynode();		//display primarynode map
	display_primarypath();		//display primarypath map
	display_pslotposition();	//display slot positions of primary lightpath
	display_backupnode();		//display backupnode map
	display_backuppath();		//display backuppath map
	display_bslotposition();	//display slot positions of backup lightpath
	display_shared_id();	//display sharing details of shared_id
	display_mainmap();		//display the whole lightpath and their slots which are in mainmap map

	display_deleted_lightpaths();	//displaying deleted lightpaths from a vector

	displaying_layer_id();		//displaying the whole id for mainnodemap and their lightpaths for each layer
	//display_mainnodemap();		//display the whole nodes and their layers which are in mainnodemap map
	display_node_lps();		//Displaying nodes and their current lightpath

	displaying_bplink_lps();	//displaying bplink and lightpaths
	display_splitmulmap();	//display link->lightpaths->splitting_lps details

		//} // if(Lightpath>=10)

} //if((No_Of_Slots!=0)&&(start!=0)&&(end!=0)&&(start!=end))
//===========================FOR REMOVING A LIGHTPATH===================================

bool lpchecker(VecMap& mypbmaps)
{
	for (VecMap::reverse_iterator revpos = mypbmaps.rbegin(); revpos != mypbmaps.rend(); ++revpos)
	{
		if (revpos->first == Lightpath)
		{
			return true;
		}
	}
	return false;
}

void printcase1()
{
	cout << "\n\nBLOCK!!!! NO PRIMARY PATH AVAILABLE (CASE-1)>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>";
}

void printcase2()
{
	cout << "\n\nBLOCK!!!! NO SLOTS AVAILABLE FOR PRIMARY PATH (CASE-2)[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[";
}

void printcase3()
{
	cout << "\n\nBLOCK!!!! NO BACKUP PATH AVAILABLE (CASE-3)$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n";
}

void printcase4()
{
	cout << "\n\nBLOCK!!!! NO SLOTS AVAILABLE FOR BACKUP PATH (CASE-4)##################################################################################";
}

void printcase5()
{
	cout << "\n\nBLOCK!!!! NO LAYERS AVAILABLE (CASE-5)---------------------------------------------------------------------\n";
}

void printcase6()
{
	cout << "\n\nBLOCK!!!! NO ADDITIONAL CONDITION SATISFIED (CASE-6)---------------------------------------------------------------------\n";
}
//=============================FOR TRIAL ENDS=========================================

//===================================IMPLEMENTING FLEX-WSS STARTS===========================================

//-----------------------------------GENERATING MAINNODEMAP STARTS----------------------------------------------
//INSERT NODES AND LAYERS INITIALLY
void insert_nodes()
{
	for (int i = 1; i <= max_nodes; i++)
	{
		mynode.insert(make_pair((i), new int[MAX_LAYERS]()));
	}
}

//COPYING LAYERS FROM MYNODE TO MAINNODEMAP
void copy_layers_mynode_mainnodemap()
{
	for (itn = mynode.begin(); itn != mynode.end(); ++itn)
	{
		for (size_t i = 0; i < MAX_LAYERS; ++i)
		{
			mainnodemap[itn->first].push_back(itn->second[i]);
		}
	}
	mynode.erase(mynode.begin(), mynode.end());
}

//DISPLAY MAINNODEMAP MAP
void display_mainnodemap()
{
	cout << "\n\n\nMAINNODEMAP CONTAINS (NODE=>LAYERS)\n";
	cout << "------------------------------";
	for (mnodemap = mainnodemap.begin(); mnodemap != mainnodemap.end(); ++mnodemap)
	{
		cout << "\n" << (*mnodemap).first << "=>";
		for (size_t n = 0; n < (*mnodemap).second.size(); n++)
		{
			cout << (*mnodemap).second[n] << "  ";
		}
	}
	cout << "\n";
}
//-----------------------------------GENERATING MAINNODEMAP ENDS--------------------------------------------------

//---------------------CHECKING LAYER FOR EACH NODE---------------------------------------------
map<int, int> ProcessLayer(int Sslot, int Eslot)
{
	IntVec BeforeVec;
	IntVec AfterVec;
	NodeLayer.clear();

	for (size_t i = 0; i < bplinks.size(), i < bpnodes.size(); i++)
	{
		//-----------------------------FOR FIRST NODE----------------------------------------
		if (bpnodes[i] == bpnodes[0])
		{
			//cout<<"\n\nNODE:"<<bpnodes[i]<<"=>"<<"LINKS:"<<bplinks[i];
			AfterVec = BeforeAfterLinks(sorted_map, bplinks[i], Sslot, Eslot);				//Getting sharable lightpaths in after node

			if (AfterVec.empty())
			{
				NodeLayer = BothNullLinkOperation(bpnodes[i]);
			}
			else if (!AfterVec.empty())
			{
				NodeLayer = AnyOneLinkOperation(AfterVec, bpnodes[i]);
			}
		}
		//-----------------------------FOR FIRST NODE----------------------------------------

		//-----------------------------FOR LAST NODE----------------------------------------
		if ((bpnodes[i] == bpnodes[bpnodes.size() - 1]))
		{
			//cout<<"\n\nNODE:"<<bpnodes[i]<<"=>"<<"LINKS:"<<bplinks[i-1];
			BeforeVec = BeforeAfterLinks(sorted_map, bplinks[i - 1], Sslot, Eslot);

			if (BeforeVec.empty())
			{
				NodeLayer = BothNullLinkOperation(bpnodes[i]);
			}
			else if (!BeforeVec.empty())
			{
				NodeLayer = AnyOneLinkOperation(BeforeVec, bpnodes[i]);
			}
		}
		//-----------------------------FOR LAST NODE----------------------------------------

		//-----------------------------FOR INTERMEDIATE NODES----------------------------------------
		if ((bpnodes[i] != bpnodes[0]) && (bpnodes[i] != bpnodes[bpnodes.size() - 1]))
		{
			BeforeVec = BeforeAfterLinks(sorted_map, bplinks[i - 1], Sslot, Eslot);			//Getting sharable lightpaths in before node

			AfterVec = BeforeAfterLinks(sorted_map, bplinks[i], Sslot, Eslot);			//Getting sharable lightpaths in after node

			if ((BeforeVec.empty()) && (AfterVec.empty()))
			{
				NodeLayer = BothNullLinkOperation(bpnodes[i]);
			}

			else if ((!BeforeVec.empty()) && (AfterVec.empty()))
			{
				NodeLayer = AnyOneLinkOperation(BeforeVec, bpnodes[i]);
			}

			else if ((BeforeVec.empty()) && (!AfterVec.empty()))
			{
				NodeLayer = AnyOneLinkOperation(AfterVec, bpnodes[i]);
			}

			else if ((!BeforeVec.empty()) && (!AfterVec.empty()))
			{
				NodeLayer = BothFullLinkOperation(BeforeVec, AfterVec, bpnodes[i]);
			}
		}
		//-----------------------------FOR INTERMEDIATE NODES----------------------------------------
	}
	return NodeLayer;
}

//-----------------BOTH ADJACENT LINKS FOR A NODE HAVE NO SHARABLE LIGHTPATHS--------------------------------------
map<int, int> BothNullLinkOperation(int node)
{
	VecMap::iterator itmain = mainnodemap.find(node); // find key equal to val
	if (itmain != mainnodemap.end())
	{
		int mainKey = itmain->first;
		IntVec& mainVec = itmain->second;
		for (size_t j = 0; j < mainVec.size(); j++)
		{
			if (mainVec[j] == 0)
			{
				NodeLayer.insert(make_pair(node, static_cast<int>(j + 1)));		//(i+1) - START LAYER FROM 1
				break;
			}
		}
	}
	return NodeLayer;
}

//-----------------ANY ONE OF THE ADJACENT LINK FOR A NODE HAS SHARABLE LIGHTPATHS--------------------------------------
map<int, int> AnyOneLinkOperation(IntVec& AnyOneLink, int node)
{
	int flag = 0;
	for (IntVec::iterator it = AnyOneLink.begin(); it != AnyOneLink.end(); ++it)
	{
		VecMap::iterator itmain = mainnodemap.find(node); // find key equal to val
		if (itmain != mainnodemap.end())
		{
			int mainKey = itmain->first;
			IntVec& mainVec = itmain->second;

			for (size_t j = 0; j < mainVec.size(); j++)
			{
				VecMap::iterator itlid = layer_id.find(mainVec[j]); // find key equal to val
				if (itlid != layer_id.end())
				{
					IntVec& LidVal = itlid->second;
					if ((IsValueInVec((*it), LidVal)))
					{
						NodeLayer.insert(make_pair(node, static_cast<int>(j + 1)));
						flag = 1;
						break;
					}
				}
			}
		}
		if (flag == 1)
		{
			break;
		}
	}

	/*
	//--------------DELETING SELECTED VECTOR FROM THE END AND CHECK FOR THE SHARING POSSIBILITIES---------------
	VecMap::iterator itmain = mainnodemap.find(node); // find key equal to val
	if (itmain != mainnodemap.end())
	{
		int mainKey=itmain->first;
		IntVec& mainVec=itmain->second;
		int count=0;
		for(size_t j=0; j<mainVec.size(); j++)
		{
			VecMap::iterator itlid = layer_id.find(mainVec[j]); // find key equal to val
			if (itlid != layer_id.end())
			{
				IntVec& LidVal=itlid->second;

				if(isVectorMatch(LidVal, AnyOneLink))
				{
					NodeLayer.insert(make_pair(node,j+1));
					break;
				}
			}
			else
			{
				count++;
			}
		}

		if(count==mainVec.size())
		{
			int flag=0;
			while(flag!=1)
			{
				if(!AnyOneLink.empty())
				{
					AnyOneLink.erase(remove(AnyOneLink.begin(), AnyOneLink.end(), (*--AnyOneLink.end())), AnyOneLink.end());
				}

				int count1=0;
				for(size_t j=0; j<mainVec.size(); j++)
				{
					VecMap::iterator itlid = layer_id.find(mainVec[j]); // find key equal to val
					if (itlid != layer_id.end())
					{
						int lidK=itlid->first;
						IntVec& LidVal=itlid->second;

						if(isVectorMatch(AnyOneLink, LidVal))
						{
							NodeLayer.insert(make_pair(node,j+1));
							flag=1;
							break;
						}
					}
					else
					{
						count1++;
					}
				}
				if(count1==mainVec.size())
				{
					flag=1;
				}
			}
		}
	}
	//--------------DELETING SELECTED VECTOR FROM THE END AND CHECK FOR THE SHARING POSSIBILITIES---------------
	*/

	return NodeLayer;
}

//--------------------BOTH ADJACENT LINKS FOR A NODE HAVE SHARABLE LIGHTPATHS---------------------

map<int, int> BothFullLinkOperation(IntVec& BeforeVec, IntVec& AfterVec, int node)
{
	int flag = 0;
	for (IntVec::iterator it1 = BeforeVec.begin(); it1 != BeforeVec.end(); ++it1)
	{
		for (IntVec::iterator it2 = AfterVec.begin(); it2 != AfterVec.end(); ++it2)
		{
			VecMap::iterator itmain = mainnodemap.find(node); // find key equal to val
			if (itmain != mainnodemap.end())
			{
				int mainKey = itmain->first;
				IntVec& mainVec = itmain->second;

				for (size_t j = 0; j < mainVec.size(); j++)
				{
					VecMap::iterator itlid = layer_id.find(mainVec[j]); // find key equal to val
					if (itlid != layer_id.end())
					{
						IntVec& LidVal = itlid->second;
						if ((IsValueInVec((*it1), LidVal)) && (IsValueInVec((*it2), LidVal)))
						{
							NodeLayer.insert(make_pair(node, static_cast<int>(j + 1)));
							flag = 1;
							break;
						}
					}
				}
			}
			if (flag == 1)
			{
				break;
			}
		}
		if (flag == 1)
		{
			break;
		}
	}
	return NodeLayer;
}

//-----------------------------------------------------------------------------

bool IsValueInVec(int value, IntVec& v2)
{
	IntVec::iterator itv = find(v2.begin(), v2.end(), value);
	if (itv != v2.end())
	{
		return true;
	}
	return false;
}

//-----------------------------------------------------------------------------

bool isVectorMatch(IntVec& v1, IntVec& v2)
{
	for (IntVec::iterator itv = v1.begin(); itv != v1.end(); ++itv)
	{
		if (!IsValueInVec((*itv), v2))
		{
			return false;
		}
	}
	return true;
}

//--------------------GENERATING SHARABLE LIGHTPATHS FROM SORTED_MAP------------------------------

IntVec BeforeAfterLinks(VecMap& sorted_map, int link, int Sslot, int Eslot)
{
	IntVec BAVec;
	BAVec.clear();

	if (!sorted_map.empty())
	{
		VecMap::iterator itsort = sorted_map.find(link);
		if (itsort != sorted_map.end())
		{
			IntVec& sortVec = itsort->second;
			//cout<<"\n\nNODE:"<<bpnodes[i]<<"=>"<<"LINKS:"<<bplinks[i-1]<<" and "<<bplinks[i];
			for (IntVec::iterator itafterVec = sortVec.begin(); itafterVec != sortVec.end(); ++itafterVec)
			{
				if (check_current_lp_overlap(Sslot, Eslot, (*itafterVec)))
				{
					BAVec.push_back(*itafterVec);
				}
			}
		}
	}
	return BAVec;
}

//--------------------DISPLAYING NODE_LAYER MAP----------------------------------------------
//DISPLAY NODE_LAYER MAP
void display_node_layer()
{
	cout << "\n\n\nNODELAYER CONTAINS (TEMPMAP) (NODE=>LAYER)\n";
	cout << "-----------------------------------------";
	for (map<int, int>::iterator it = NodeLayer.begin(); it != NodeLayer.end(); ++it)
	{
		cout << "\n" << it->first << "=>" << it->second;
	}
	cout << "\n";
}
//--------------------DISPLAYING NODE_LAYER MAP----------------------------------------------

//---------------------INSERTING LIGHTPATHS INTO LAYERS_ID, MAINNODEMAP MAPS---------------------------------------

void insert_lightpath_into_layer()
{
	for (size_t i = 0; i < bpnodes.size(); i++)
	{
		map<int, int>::iterator itNL = NodeLayer.find(bpnodes[i]);
		if (itNL != NodeLayer.end())

			int node = itNL->first;
		int layer = itNL->second;

		VecMap::iterator itmainnode = mainnodemap.find(bpnodes[i]);
		if (itmainnode != mainnodemap.end())
		{
			IntVec& mainnodeVec = itmainnode->second;

			int key = 0;
			if (mainnodeVec[layer - 1] == 0)
			{
				if (layer_id.empty())
				{
					key = LAYERID;
				}
				else
				{
					VecMap::iterator it = layer_id.end();
					it--;
					key = it->first + 1; // new key
				}
				layer_id[key].push_back(Lightpath);
				mainnodeVec[layer - 1] = key;
			}
			else
			{
				layer_id[mainnodeVec[layer - 1]].push_back(Lightpath);
			}
		}
	}
}
//---------------------INSERTING LIGHTPATHS INTO LAYERS_ID, MAINNODEMAP MAPS---------------------------------------

//---------------------DELETING LIGHTPATHS FROM LAYERS---------------------------------------
void del_lightpath_in_layers(int dLightpath)
{
	VecMap::iterator iback = backupnode.find(dLightpath);
	if (iback != backupnode.end())
	{
		IntVec& delnodes = iback->second;

		if (!delnodes.empty())		//Checks if dLightpath available in mypathbackup
		{
			for (IntVec::iterator itdel = delnodes.begin(); itdel != delnodes.end(); ++itdel)

			{
				VecMap::iterator imain = mainnodemap.find(*itdel);
				if (imain != mainnodemap.end())
				{
					IntVec& mainVec = imain->second;
					for (size_t i = 0; i < mainVec.size(); i++)
					{
						if (mainVec[i] != 0)
						{
							for (VecMap::iterator itlid = layer_id.begin(); itlid != layer_id.end(); )
							{
								IntVec& lidVec = itlid->second;

								if (itlid->first == mainVec[i])
								{

									lidVec.erase(remove(lidVec.begin(), lidVec.end(), (dLightpath)), lidVec.end());

									if (lidVec.size() == 0) // modified vector has no entry
									{
										layer_id.erase(itlid++);   // remove entry from 'sorted_map'
										mainVec[i] = 0;
									}
									else
									{
										++itlid;
									}
								}
								else
								{
									++itlid;
								}
							}
						}
					}
				}
			}
		}
	}
}
//---------------------DELETING LIGHTPATHS FROM LAYERS---------------------------------------

//------------------------DISPLAYING LAYER_ID-----------------------------------------------
void displaying_layer_id()
{
	cout << "\n\nLAYER_ID CONTAINS (LAYER-ID => LIGHTPATHS\n";
	cout << "----------------------------------------";
	for (VecMap::iterator it = layer_id.begin(); it != layer_id.end(); ++it)
	{
		cout << endl << it->first << " =>";
		for (size_t n = 0; n < (*it).second.size(); n++)
		{
			cout << " " << it->second[n];
		}
	}
	cout << "\n\n";
}
//------------------------DISPLAYING LAYER_ID-----------------------------------------------

//-------------------INSERTING LIGHTPATHS INTO NODE_LPS MAP---------------------------------------
void insert_lightpath_into_node_lps()
{
	for (size_t i = 0; i < bpnodes.size(); i++)
	{
		node_lps[bpnodes[i]].insert(Lightpath);
	}
}
//-------------------INSERTING LIGHTPATHS INTO NODE_LPS MAP---------------------------------------

//-------------------DELETING LIGHTPATHS FROM NODE_LPS MAP---------------------------------------
void del_lightpath_in_node_lps(int dLightpath)
{
	VecMap::iterator iback = backupnode.find(dLightpath);
	if (iback != backupnode.end())
	{
		IntVec& delnodes = iback->second;

		if (!delnodes.empty())		//Checks if dLightpath available in mypathbackup
		{
			for (IntVec::iterator itd = delnodes.begin(); itd != delnodes.end(); ++itd)
			{
				for (map<int, set<int>>::iterator it = node_lps.begin(); it != node_lps.end(); )
				{
					int key = it->first;

					if ((*itd) == key)
					{
						set<int>& itVal = it->second;
						itVal.erase(dLightpath);

						if (itVal.size() == 0) // modified vector has no entry
						{
							node_lps.erase(it++);   // remove entry from 'sorted_map'
						}
						else
						{
							++it;
						}
					}
					else
					{
						++it;
					}
				}
			}
		}
	}
}

//-------------------DISPLAYING NODE-LIGHTPATHS IN NODE_LPS MAP---------------------------------------
void display_node_lps()
{
	cout << "\n\n\nNODE_LIGHTPATHS CONTAINS (NODE=>LIGHTPATHS)\n";
	cout << "-----------------------------------";
	for (map<int, set<int>>::iterator it = node_lps.begin(); it != node_lps.end(); ++it)
	{
		cout << endl << it->first << " =>";
		for (set<int>::iterator sit = it->second.begin(); sit != it->second.end(); ++sit)
		{
			cout << " " << (*sit);
		}
	}
	cout << "\n\n";
}
//-------------------DISPLAYING NODE-LIGHTPATHS IN NODE_LPS MAP---------------------------------------

//========================IMPLEMENTING FLEX-WSS ENDS===========================================


//==============================IMPLEMENTING PBPS STARTS===========================================
//---------------------------------LINKS AND THEIR LIGHTPATHS------------------------------------
//-----------------------------INSERTING LINK AND LIGHTPATHS INTO BPLINK_LPS--------------------------------
void insert_lightpath_into_link_lps()
{
	for (IntVec::iterator it = bplinks.begin(); it != bplinks.end(); ++it)
	{
		bplink_lps[(*it)].insert(Lightpath);
	}
}
//-----------------------------INSERTING LINK AND LIGHTPATHS INTO LINK_LPS--------------------------------

//------------------------DISPLAYING LINK_LPS-----------------------------------------------
void displaying_bplink_lps()
{
	cout << "\n\n\nBPLINK_LPS (BACKUP LINKS => LIGHTPATHS\n";
	cout << "------------------------------------";
	for (MapSet::iterator it = bplink_lps.begin(); it != bplink_lps.end(); ++it)
	{
		cout << endl << it->first << " =>";
		for (IntSet::iterator its = (*it).second.begin(); its != (*it).second.end(); ++its)
		{
			cout << " " << (*its);
		}
	}
	cout << "\n\n";
}
//------------------------DISPLAYING LINK_LPS-----------------------------------------------

//--------------------------DELETING LIGHTPATHS FROM LINK_LPS--------------------------------
void del_lightpath_from_bplink_lps(IntVec& bplVec, int dLightpath)
{
	for (IntVec::iterator it = bplVec.begin(); it != bplVec.end(); ++it)
	{
		for (MapSet::iterator itk = bplink_lps.begin(); itk != bplink_lps.end(); )
		{
			int key = itk->first;
			IntSet& dellps = itk->second;

			if (key == (*it))
			{
				assert(!dellps.empty()); // check there is something before end
				IntSet::iterator itv = dellps.end();
				--itv;

				dellps.erase(dLightpath);  // remove 'value' from vector

				if (dellps.size() == 0) // modified vector has no entry
				{
					bplink_lps.erase(itk++);   // remove entry from 'sorted_map'
				}
				else
				{
					++itk;
				}
			}
			else
			{
				++itk;
			}
		}
	}
}
//--------------------------DELETING LIGHTPATHS FROM LINK_LPS--------------------------------

/*
//--------------------------DELETING LIGHTPATHS FROM LINK_LPS (IN TRIAL USAGE)--------------------------------
void del_lightpath_from_bplink_lps_1(int Lightpath)
{
	for(MapSet::iterator itk = bplink_lps.begin(); itk != bplink_lps.end(); )
	{
		int key=itk->first;
		IntSet & dellps=itk->second;

		assert(!dellps.empty()); // check there is something before end
		IntSet::iterator itv = dellps.end();
		--itv;

		dellps.erase(remove(dellps.begin(), dellps.end(), Lightpath), dellps.end());  // remove 'value' from vector

		if (dellps.size() == 0) // modified vector has no entry
		{
			bplink_lps.erase(itk++);   // remove entry from 'sorted_map'
		}
		else
		{++itk;}

	}
}
//--------------------------DELETING LIGHTPATHS FROM LINK_LPS (IN TRIAL USAGE)--------------------------------
*/
//---------------------------------LINKS AND THEIR LIGHTPATHS------------------------------------

//-------------------------------SPLITTING DETAILS(INSERTING, DISPLAYING, DELETING)----------------------
//----------------------------------------INSERTING---------------------------------------------------
bool NotFind(set<int>& myset, int value)
{
	for (set<int>::iterator sSet = myset.begin(); sSet != myset.end(); ++sSet)
	{
		if (value == (*sSet))
		{
			return false;
		}
	}
	return true;
}

//----------------------------------------------------------------------------------------
void insert_splitmulmap(const IntVec& backvec)
{
	IntSet splitSet;		//To store the current lightpath and shared lightpath (related with splitMap)
	splitSet.clear();
	MapSet splitMap;		//To store the links => current lightpath and shared lightpath (related with splitVec)
	splitMap.clear();
	IntSet TempSet;
	bsltposn = backslotposition.find(Lightpath);
	if (bsltposn != backslotposition.end())
	{
		IntVec& posvec = bsltposn->second;
		int start = posvec[0];
		int end = posvec[1];
		for (itmm = mainmap.begin(); itmm != mainmap.end(); ++itmm)
		{
			int itmmkey = itmm->first;
			IntVec& itmmvec = itmm->second;
			IntVec::const_iterator itv = find(backvec.begin(), backvec.end(), itmmkey);
			if (itv != backvec.end())
			{
				for (int index = start; index <= end; index++)
				{
					VecMap::iterator its = shared_id.find(itmmvec[index]); // find key equal to val
					if (its != shared_id.end())
					{
						IntVec& cVec = its->second; // note reference

						copy(cVec.begin(), cVec.end(), inserter(splitSet, splitSet.end()));	//Copy a vector to a set

						for (set<int>::iterator sptSet = splitSet.begin(); sptSet != splitSet.end(); ++sptSet)
							if ((*sptSet) != Lightpath)
							{
								splitMap[Lightpath].insert(*sptSet);
							}
					}
				}

				TempSet.clear();
				for (int index = start; index <= end; index++)
				{
					VecMap::iterator its = shared_id.find(itmmvec[index]); // find key equal to val
					if (its != shared_id.end())
					{
						if (NotFind(TempSet, itmmkey))
						{
							SplitMulMap.insert(make_pair(itmmkey, splitMap));		//for(size_t n=0; n<bplinks.size(); n++) {if(myKey==bplinks[n]) {SplitMap.insert(make_pair(myKey, splitMap))}} "for nodes instead of links"
							TempSet.insert(itmmkey);
						}
					}
				}
			}
		}
	}
}
//-------------------------------------------INSERTING-----------------------------------------------------

//------------------------------------------DISPLAYING------------------------------------------------------
void display_splitmulmap()
{
	cout << "\n\nSPLITMULMAP CONTAINS (LINK=>CUR.LP----->SPLITTING LPS\n";
	cout << "---------------------------------------------";
	for (TriMul::iterator it = SplitMulMap.begin(); it != SplitMulMap.end(); ++it)
	{
		cout << endl << it->first << " =>";
		MapSet& Sec = it->second;
		for (MapSet::iterator itm = Sec.begin(); itm != Sec.end(); ++itm)
		{
			cout << itm->first << " ---->";
			for (IntSet::iterator itmm = (*itm).second.begin(); itmm != (*itm).second.end(); ++itmm)
			{
				cout << " " << (*itmm);
			}
		}
	}
	cout << "\n\n";
}
//------------------------------------------DISPLAYING------------------------------------------------------

//------------------------------------------DELETING------------------------------------------------------
void del_lightpath_from_splitmulmap(int dLightpath)
{
	for (TriMul::iterator it = SplitMulMap.begin(); it != SplitMulMap.end(); ++it)
	{
		MapSet& Sec = it->second;

		for (MapSet::iterator itm = Sec.begin(); itm != Sec.end(); )
		{
			int key = itm->first;
			IntSet& SetVal = itm->second;

			assert(!SetVal.empty()); // check there is something before end
			IntSet::iterator its = SetVal.end();
			--its;

			if (key == dLightpath)
			{
				Sec.erase(itm++);
			}
			else
			{
				SetVal.erase(dLightpath);
				if (SetVal.size() == 0) // modified vector has no entry
				{
					Sec.erase(itm++);   // remove entry from 'sorted_map'
				}
				else
				{
					{++itm; }
				}
			}
		}
	}

	for (TriMul::iterator it = SplitMulMap.begin(); it != SplitMulMap.end();)
	{
		MapSet& mySec = it->second;
		if (mySec.size() == 0) // modified vector has no entry
		{
			SplitMulMap.erase(it++);   // remove entry from 'sorted_map'
		}
		else
		{
			{++it; }
		}
	}
}
//------------------------------------------DELETING------------------------------------------------------
//-------------------------------SPLITTING DETAILS(INSERTING, DISPLAYING, DELETING)----------------------

//=================================PBPS CONDITION-1 MAIN FUNCTION==========================================

IntSet PBPS_COND_1(int start, int end)
{
	NotifySet1.clear();
	{

		for (size_t i = 0; i < bplinks.size(); i++)
		{
			if (bplinks[i] != bplinks[0])
			{
				if (FindSet(FourLinks, bplinks[i]))
				{
					if (bplink_lps.empty())
					{
						copy(bplinks.begin(), bplinks.end(), inserter(NotifySet1, NotifySet1.end()));	//Copy a vector to a set
						NotifySet1.erase(NotifySet1.begin());
						//NotifySet1.insert(bplinks[i]);
					}
					else
					{
						MapSet::iterator ibp2 = bplink_lps.find(bplinks[i]);
						if (ibp2 != bplink_lps.end())
						{
							IntSet& bplpVec2 = ibp2->second;
							//assert(!bplpVec2.empty()); // check there is something before end
							//IntSet::iterator itv2 = bplpVec2.end();
							//--itv2;

							if (SplitMulMap.empty())
							{
								copy(bplinks.begin(), bplinks.end(), inserter(NotifySet1, NotifySet1.end()));	//Copy a vector to a set
								NotifySet1.erase(NotifySet1.begin());
							}
							else
							{
								for (IntSet::iterator itv2 = bplpVec2.begin(); itv2 != bplpVec2.end(); ++itv2)
								{
									TriMul::iterator itri = SplitMulMap.find(bplinks[i]);
									if (itri != SplitMulMap.end())
									{
										int itmulkey = itri->first;
										MapSet& SplFirst = itri->second;

										MapSet::iterator itmul = SplFirst.find((*itv2));
										if (itmul != SplFirst.end())
										{
											IntSet& SplThird = itmul->second;

											for (IntSet::iterator itTVec = SplThird.begin(); itTVec != SplThird.end(); ++itTVec)
											{
												if (!check_current_lp_overlap(start, end, (*itTVec)))
												{
													NotifySet1.insert(bplinks[i]);
												}
												else //if(check_current_lp_overlap(start, end, (*itTVec)))
												{
													if (!check_is_curr_joint(bpnodes, (*itv2)))
													{
														NotifySet1.insert(bplinks[i]);
													}
												}
											}
										}
										else		//MapSet::iterator itmul = SplFirst.find((*itv2)); 
										{
											NotifySet1.insert(bplinks[i]);
										}
									}
									else		//TriMul::iterator itri = SplitMulMap.find(bplinks[i]); 
									{
										NotifySet1.insert(bplinks[i]);
									}
								}
							}
						}
						else		//MapSet::iterator ibp2 = bplink_lps.find(bplinks[i]); 
						{
							NotifySet1.insert(bplinks[i]);
						}
					}
				}
				else		//MapSet::iterator ibp2 = bplink_lps.find(bplinks[i]); 
				{
					NotifySet1.insert(bplinks[i]);
				}
			}
		}
	}
	return NotifySet1;
}

//-------------------------------------------------------------------------------------
// returns true if finds value in 'set'
bool FindSet(const IntSet& iSet, int value)
{
	return find(iSet.begin(), iSet.end(), value) != iSet.end();
}

//=================================PBPS CONDITION-2 MAIN FUNCTION==========================================
IntSet PBPS_COND_2(int start, int end)
{
	NotifySet2.clear();
	{

		for (size_t i = 0; i < bplinks.size(); i++)
		{
			if (bplinks[i] != bplinks[0])
			{
				if (FindSet(FourLinks, bplinks[i]))
				{
					if (bplink_lps.empty())
					{
						copy(bplinks.begin(), bplinks.end(), inserter(NotifySet2, NotifySet2.end()));	//Copy a vector to a set
						NotifySet2.erase(NotifySet2.begin());
						//NotifySet2.insert(bplinks[i]);
					}
					else
					{

						MapSet::iterator ibp2 = bplink_lps.find(bplinks[i]);
						if (ibp2 != bplink_lps.end())
						{
							IntSet& bplpVec2 = ibp2->second;
							assert(!bplpVec2.empty()); // check there is something before end
							IntSet::iterator itv2 = bplpVec2.end();
							--itv2;

							//for(IntSet::iterator itv2 = bplpVec2.begin(); itv2 != bplpVec2.end(); ++itv2)

							MapSet::iterator ibp1 = bplink_lps.find(bplinks[i - 1]);
							if (ibp1 != bplink_lps.end())
							{
								IntSet& bplpVec1 = ibp1->second;
								assert(!bplpVec1.empty()); // check there is something before end
								IntSet::iterator itv1 = bplpVec1.end();
								--itv1;

								if (!FindSet(bplpVec2, (*itv1)))
								{
									VecMap::iterator itsort = sorted_map.find(bplinks[i - 1]);
									if (itsort != sorted_map.end())
									{
										IntVec& sortvec = itsort->second;
										for (IntVec::iterator itS = sortvec.begin(); itS != sortvec.end(); ++itS)
										{
											if (check_current_lp_overlap(start, end, (*itS)))			//Current and previous node lightpath overlap
											{
												if ((!check_overlap_bp_slots((*itv1), (*itv2))))
												{
													{NotifySet2.insert(bplinks[i]); }
												}
												else //if((check_overlap_bp_slots((*itv1), (*itv2))))
												{
													if (!check_is_joint((*itv1), (*itv2)))
													{
														{NotifySet2.insert(bplinks[i]); }
													}
												}
											}
											else		//if(check_current_lp_overlap(start, end, (*itS)))
											{
												NotifySet2.insert(bplinks[i]);
											}
										}
									}
									else		//VecMap::iterator itsort = sorted_map.find(*itv1); 
									{
										NotifySet2.insert(bplinks[i]);
									}
								}
								else		//if((!Find(bplpVec2, (*itv1))) 
								{
									NotifySet2.insert(bplinks[i]);
								}
							}
							else		//MapSet::iterator ibp1 = bplink_lps.find(bplinks[i-1]);  
							{
								NotifySet2.insert(bplinks[i]);
							}
						}
						else		//MapSet::iterator ibp2 = bplink_lps.find(bplinks[i]); 
						{
							NotifySet2.insert(bplinks[i]);
						}
					}
				}
				else		//MapSet::iterator ibp2 = bplink_lps.find(bplinks[i]); 
				{
					NotifySet2.insert(bplinks[i]);
				}
			}
		}
	}
	return NotifySet2;
}

/*
--------------------------CHECKING PBPS'S BOTH CONDITIONS SIMULTANEOUSLY (ERROR IN THIS CODE)----------
set<int> PBPS(int start, int end)
{
NotifySet.clear();
//if(Lightpath!=0)
{
	//if(bplink_lps.empty())
	//{
		//NotifySet.insert(1);
	//}
	//else
	for(size_t i=0; i<bplinks.size(); i++)
	{
	if(bplinks[i]!=bplinks[0])
	{
		VecMap::iterator ibp2 = bplink_lps.find(bplinks[i]);
		if (ibp2 != bplink_lps.end())
		{
			IntVec& bplpVec2=ibp2->second;
			assert(!bplpVec2.empty()); // check there is something before end
			IntVec::iterator itv2 = bplpVec2.end();
			--itv2;

			if(!bplpVec2.empty())
			{
			VecMap::iterator ibp1 = bplink_lps.find(bplinks[i-1]);
			if (ibp1 != bplink_lps.end())
			{
				IntVec& bplpVec1=ibp1->second;
				assert(!bplpVec1.empty()); // check there is something before end
				IntVec::iterator itv1 = bplpVec1.end();
				--itv1;

				TriMul::iterator itri = SplitMulMap.find(bplinks[i]);
				if (itri != SplitMulMap.end())
				{
					int itmulkey=itri->first;
					MapSet& SplFirst=itri->second;

					MapSet::iterator itmul = SplFirst.find((*itv2));
					if (itmul != SplFirst.end())
					{
						IntSet& SplThird=itmul->second;

						for(IntSet::iterator itTVec = SplThird.begin(); itTVec != SplThird.end(); ++itTVec)
						{
							if(check_current_lp_overlap(start, end, (*itTVec)))
							{
								NotifySet.insert(2); NotifySet.insert(2);
							}
						}
					}
				}

				if((Find(bplpVec1, (*itv2)))&&(check_overlap_bp_slots((*itv1), (*itv2))))
				{
					if(check_disjoint ((*itv1), (*itv2)))
					{
						NotifySet.insert(3); NotifySet.insert(3); NotifySet.insert(3);
					}
				}
			}
			}
		}
		//else
		//{NotifySet.insert(1);}
	}
}
}
return NotifySet;
}
*/

//========================For checking current lp with other splitting history (CASE-1)===========================================
bool check_is_curr_joint(IntVec& curVec, int lp)
{
	VecMap::iterator itp1 = primarynode.find(lp);
	if (itp1 != primarynode.end())
	{
		IntVec& bpVec = itp1->second;

		if ((vmatch(curVec, bpVec)))
		{
			return true;
		}
	}
	return false;
}

//-----------------------------------------------------------------------------------------------

bool check_current_lp_overlap(int Sslot, int Eslot, int lp)
{
	VecMap::iterator it = backslotposition.find(lp);
	if (it != backslotposition.end())
	{
		IntVec& vec = it->second;

		for (int i = Sslot; i <= Eslot; i++)
		{
			for (int j = vec[0]; j <= vec[1]; j++)
			{
				if (i == j)
				{
					return true;
				}
			}
		}
	}
	return false;
}

//========================For checking current lp with other splitting history (CASE-1)===========================================

//====================For checking other lps with current lp's splitting history (CASE-2)=========================================

bool check_is_joint(int lp1, int lp2)
{
	VecMap::iterator itp1 = primarynode.find(lp1);
	if (itp1 != primarynode.end())
	{
		IntVec& curVec = itp1->second;

		VecMap::iterator itp2 = primarynode.find(lp2);
		if (itp2 != primarynode.end())
		{
			IntVec& bpVec = itp2->second;

			if ((vmatch(curVec, bpVec)))
			{
				return true;
			}
		}
	}
	return false;
}

//-----------------------------------------------------------------------------------------------

bool check_overlap_bp_slots(int lp1, int lp2)
{
	VecMap::iterator it1 = backslotposition.find(lp1);
	if (it1 != backslotposition.end())
	{
		IntVec& vec1 = it1->second;

		VecMap::iterator it2 = backslotposition.find(lp2);
		if (it2 != backslotposition.end())
		{
			IntVec& vec2 = it2->second;

			for (int i = vec1[0]; i <= vec1[1]; i++)
			{
				for (int j = vec2[0]; j <= vec2[1]; j++)
				{
					if (i == j)
					{
						return true;
					}
				}
			}
		}
	}
	return false;
}
//====================For checking other lps with current lp's splitting history (CASE-2)===================================

//===================================IMPLEMENTING PBPS ENDS===========================================

//=============================FOR AUTOMATION STARTS====================================

//=====following random number-generators are of UNIFORM DISTRIBUTION==============
//will create including 'lowest_number' and 'highest_number'.

int random_source(int lowest_number, int highest_number)
{
	int range = highest_number - lowest_number + 1;
	//cout<<"\n lowest number:"<<lowest_number<<"\t"<<"highest\t"<<highest_number<<"range\t"<<range<<"RAND_MAX\t"<<RAND_MAX;
	//cout<<"\n\t\t**********"<<(lowest_number + int((float)range * (float)rand()/(RAND_MAX + 1.0)));
	return (lowest_number + int((float)range * (float)rand() / (RAND_MAX + 1.0)));
}

int random_destination(int lowest_number, int highest_number, int Source)//creates not-Source!
{
	int range = highest_number - lowest_number + 1;
	int Random_No;

	do {
		Random_No = lowest_number + int((float)range * (float)rand() / (RAND_MAX + 1.0));
	} while (Random_No == Source);

	return Random_No;
}

//=====following random number-generators are of EXPONENTIAL DISTRIBUTION==============
double expon(double mean)			//returns random exponential-variate with specified mean
{
	double u;
	do
	{
		u = rand() / (RAND_MAX + 1.0); //values in between 0 and 1.
	} while (u == 0);
	return(-mean * log(u));
}

//DISPLAY TIME SCALE AND LIGHTPATH FROM TIME_SCALE MULTIMAP		
void display_TIME_SCALE()
{
	cout << "\n\nTIME_SCALE CONTAINS (TIME=>LIGHTPATH)\n";
	cout << "---------------------------------";
	for (multimap<double, int>::iterator it = TIME_SCALE.begin(); it != TIME_SCALE.end(); ++it)
	{
		cout << endl << it->first << " => " << it->second;
	}
	cout << "\n";
}
//=============================FOR AUTOMATION ENDS====================================

//=============================CALCULATION==============================================
void Calculation()
{

	//=========================FOR TRACKING CURRENT TIME STARTS===========================
	time(&tim);
	char mytime[123];
	ctime_s(mytime, _countof(mytime), &tim);
	//=========================FOR TRACKING CURRENT TIME STARTS===========================

	double bandblock;
	double bprob;
	int totslotsused = 0;
	double spec_effe;

	//int XBlock;

	cout << "\n\n==============================================================";
	cout << "\n---------------------SIMULATION RESULTS-----------------------";
	cout << "\n-----------------------(PBPS Results)-------------------------";
	cout << "\n==============================================================";

	cout << "\n\nNetwork Graph: Bi-directinal Links-" << max_edges / 2 << "; Nodes-" << max_nodes;

	cout << "\n\nData rates used: ";
	for (multimap<int, int>::iterator itdr = DataRatemm.begin(); itdr != DataRatemm.end(); ++itdr)
	{
		cout << (itdr->first) << "Gbps" << " ";
	}

	cout << "\n\nLoad: " << LOAD << " Erlang";
	cout << "\n\nMean Holding Time: " << MeanHoldTime;
	cout << "\n\nMean Arrival Time: " << MeanArrivalTime;

	cout << "\n\nTotal number of Frequency slots: " << MAX_SLOTS;
	cout << "\n\nTotal number of Lightpaths: " << MAX_LIGHTPATHS;
	cout << "\n\nTotal number of Lightpaths departed: " << DelVec.size();
	cout << "\n\nTotal number of Layers: " << MAX_LAYERS;

	//cout<<"\n\n(Actual total load (in Gbps) is: "<<total_load()<<")";

	//cout<<"\n\nTotal load (in Gbps): "<<TOT_LOAD;

	cout << "\n\nTotal number of Lightpaths blocked: " << MAX_LIGHTPATHS - (backupnode.size() + DelVec.size());

	bprob = (double(MAX_LIGHTPATHS) - (double((backupnode.size()) + double(DelVec.size())))) / double(MAX_LIGHTPATHS);
	cout << "\n\nBlocking Probability: " << bprob;

	cout << "\n\nBlocking Details";
	cout << "\n-----------------";
	cout << "\n1-No Primary Path: " << NOPRIPATH.size();
	cout << "\n\n2-No Primary Slots: " << NOPRISLOT.size();
	cout << "\n\n3-No Backup Path: " << NOBACKPATH.size();
	cout << "\n\n4-No Backup Slots: " << NOBACKSLOT.size();
	cout << "\n\n5-No Layers: " << NOLAYER.size();
	cout << "\n\n6-No Additional Condition Satisfied: " << NOADDCON.size();

	//XBlock=(MAX_LIGHTPATHS-(backupnode.size()+DelVec.size()))-(NOPRIPATH.size()+NOPRISLOT.size()+NOBACKPATH.size()+NOBACKSLOT.size()+NOLAYER.size()+NOADDCON.size());
	cout << "\n\n7-Either No Backup Slots/No Layers/No Additional Condition or no backup path (For the trials): " << XBLOCK.size();
	cout << "\n\n(Total Blocking): " << NOPRIPATH.size() + NOPRISLOT.size() + NOBACKPATH.size() + NOBACKSLOT.size() + NOLAYER.size() + NOADDCON.size() + XBLOCK.size();

	bandblock = ((double(count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 10)) * 10) + (double(count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 40)) * 40) + (double(count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 100)) * 100) + (double(count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 400)) * 400) + (double(count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 1000)) * 1000)) / (double(RATETOTAL[0] * 10) + double(RATETOTAL[1] * 40) + double(RATETOTAL[2] * 100) + double(RATETOTAL[3] * 400) + double(RATETOTAL[4] * 1000));
	cout << "\n\nBandwidth Blocking Probability: " << bandblock;

	cout << "\n\nBlocked(10): " << (count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 10));
	cout << "\nBlocked(40): " << (count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 40));
	cout << "\nBlocked(100): " << (count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 100));
	cout << "\nBlocked(400): " << (count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 400));
	cout << "\nBlocked(1000): " << (count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 1000));

	cout << "\n\nTotal(10): " << RATETOTAL[0];
	cout << "\nTotal(40): " << RATETOTAL[1];
	cout << "\nTotal(100): " << RATETOTAL[2];
	cout << "\nTotal(400): " << RATETOTAL[3];
	cout << "\nTotal(1000): " << RATETOTAL[4];

	cout << "\n\nNo. of Slots used";
	cout << "\n------------------\n";
	for (VecMap::iterator itvm = mainmap.begin(); itvm != mainmap.end(); ++itvm)
	{
		IntVec& countVecm = itvm->second;
		cout << itvm->first << "->" << MAX_SLOTS - (count(countVecm.begin(), countVecm.end(), 0)) << "\n";
		totslotsused = totslotsused + MAX_SLOTS - static_cast<int>(count(countVecm.begin(), countVecm.end(), 0));
	}

	cout << "\nTotal slots used in the network: " << totslotsused << "\n";

	cout << "\n\nTotal spectrum used in the network: " << totslotsused * FIXED << " GHz" << "\n";	//*50 for fixed grid

	spec_effe = (double)(RATETOTAL[0] * 10 + RATETOTAL[1] * 40 + RATETOTAL[2] * 100 + RATETOTAL[3] * 400 + RATETOTAL[4] * 1000) / (double)(totslotsused * FIXED);

	cout << "\n\nSpectrum efficiency at " << (RATETOTAL[0] * 10 + RATETOTAL[1] * 40 + RATETOTAL[2] * 100 + RATETOTAL[3] * 400 + RATETOTAL[4] * 1000) << " Gbps is " << spec_effe << " bit/s/Hz";

	cout << "\n\nNo. of Layers used";
	cout << "\n------------------\n";
	for (VecMap::iterator itvec = mainnodemap.begin(); itvec != mainnodemap.end(); ++itvec)
	{
		IntVec& countVec = itvec->second;
		cout << itvec->first << "->" << MAX_LAYERS - (count(countVec.begin(), countVec.end(), 0)) << "\n";
	}

	cout << "\n\n\nSimulation End Time: " << mytime << "\n";

	cout << "\n\n----------------------END OF SIMULATION-----------------------\n";

	totslotsused = 0;

	fout << "\n==============================================================";
	fout << "\n---------------------SIMULATION RESULTS-----------------------";
	fout << "\n-----------------------(PBPS Results)-------------------------";
	fout << "\n==============================================================";

	fout << "\nNetwork Graph: Bi-directinal Links-" << max_edges / 2 << "; Nodes-" << max_nodes;

	fout << "\nData rates used: ";
	for (multimap<int, int>::iterator itdr = DataRatemm.begin(); itdr != DataRatemm.end(); ++itdr)
	{
		fout << (itdr->first) << "Gbps" << " ";
	}

	fout << "\nLoad: " << LOAD << " Erlang";
	fout << "\nMean Holding Time: " << MeanHoldTime;
	fout << "\nMean Arrival Time: " << MeanArrivalTime;

	fout << "\nTotal number of Frequency slots: " << MAX_SLOTS;
	fout << "\nTotal number of Lightpaths: " << MAX_LIGHTPATHS;

	//fout<<"\n\n(Actual total load (in Gbps) is: "<<total_load()<<")";

	//fout<<"\n\nTotal load (in Gbps): "<<TOT_LOAD;

	fout << "\nTotal number of Lightpaths departed: " << DelVec.size();
	fout << "\nTotal number of Layers: " << MAX_LAYERS;
	fout << "\nTotal number of Lightpaths blocked: " << MAX_LIGHTPATHS - (backupnode.size() + DelVec.size());

	bprob = (double(MAX_LIGHTPATHS) - (double((backupnode.size()) + double(DelVec.size())))) / double(MAX_LIGHTPATHS);
	fout << "\nBlocking Probability: " << bprob;

	fout << "\n\nBlocking Details";
	fout << "\n-----------------";
	fout << "\n1-No Primary Path: " << NOPRIPATH.size();
	fout << "\n2-No Primary Slots: " << NOPRISLOT.size();
	fout << "\n3-No Backup Path: " << NOBACKPATH.size();
	fout << "\n4-No Backup Slots: " << NOBACKSLOT.size();
	fout << "\n5-No Layers: " << NOLAYER.size();
	fout << "\n6-No Additional Condition Satisfied: " << NOADDCON.size();

	//XBlock=(MAX_LIGHTPATHS-(backupnode.size()+DelVec.size()))-(NOPRIPATH.size()+NOPRISLOT.size()+NOBACKPATH.size()+NOBACKSLOT.size()+NOLAYER.size()+NOADDCON.size());
	fout << "\n7-Either No Backup Slots/No Layers/No Additional Condition or no backup path (For the trials): " << XBLOCK.size();
	fout << "\n(Total Blocking): " << NOPRIPATH.size() + NOPRISLOT.size() + NOBACKPATH.size() + NOBACKSLOT.size() + NOLAYER.size() + NOADDCON.size() + XBLOCK.size();

	bandblock = ((double(count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 10)) * 10) + (double(count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 40)) * 40) + (double(count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 100)) * 100) + (double(count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 400)) * 400) + (double(count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 1000)) * 1000)) / (double(RATETOTAL[0] * 10) + double(RATETOTAL[1] * 40) + double(RATETOTAL[2] * 100) + double(RATETOTAL[3] * 400) + double(RATETOTAL[4] * 1000));
	fout << "\n\nBandwidth Blocking Probability: " << bandblock;

	fout << "\n\nBlocked(10): " << (count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 10));
	fout << "\nBlocked(40): " << (count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 40));
	fout << "\nBlocked(100): " << (count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 100));
	fout << "\nBlocked(400): " << (count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 400));
	fout << "\nBlocked(1000): " << (count(RATESBLOCKED.begin(), RATESBLOCKED.end(), 1000));

	fout << "\n\nTotal(10): " << RATETOTAL[0];
	fout << "\nTotal(40): " << RATETOTAL[1];
	fout << "\nTotal(100): " << RATETOTAL[2];
	fout << "\nTotal(400): " << RATETOTAL[3];
	fout << "\nTotal(1000): " << RATETOTAL[4];

	fout << "\n\nNo. of Slots used";
	fout << "\n------------------\n";
	for (VecMap::iterator itvm = mainmap.begin(); itvm != mainmap.end(); ++itvm)
	{
		IntVec& countVecm = itvm->second;
		fout << itvm->first << "->" << MAX_SLOTS - (count(countVecm.begin(), countVecm.end(), 0)) << "\n";
		totslotsused = totslotsused + MAX_SLOTS - static_cast<int>(count(countVecm.begin(), countVecm.end(), 0));
	}

	fout << "\nTotal slots used in the network: " << totslotsused << "\n";

	fout << "\n\nTotal spectrum used in the network: " << totslotsused * FIXED << " GHz" << "\n";	//*50 for fixed grid

	spec_effe = (double)(RATETOTAL[0] * 10 + RATETOTAL[1] * 40 + RATETOTAL[2] * 100 + RATETOTAL[3] * 400 + RATETOTAL[4] * 1000) / (double)(totslotsused * FIXED);

	fout << "\n\nSpectrum efficiency at " << (RATETOTAL[0] * 10 + RATETOTAL[1] * 40 + RATETOTAL[2] * 100 + RATETOTAL[3] * 400 + RATETOTAL[4] * 1000) << " Gbps is " << spec_effe << " bit/s/Hz";

	fout << "\n\nNo. of Layers used";
	fout << "\n------------------\n";
	for (VecMap::iterator itvec = mainnodemap.begin(); itvec != mainnodemap.end(); ++itvec)
	{
		IntVec& countVec = itvec->second;
		fout << itvec->first << "->" << MAX_LAYERS - (count(countVec.begin(), countVec.end(), 0)) << "\n";
	}

	fout << "\nSimulation Finish Time: " << mytime << "\n";

	/*
	for(IntSet::iterator it1=NOPRIPATH.begin(); it1!=NOPRIPATH.end(); ++it1)
	{
		fout<<(*it1)<<"  ";
	}
	fout<<"\n\n";
	for(IntSet::iterator it2=NOPRISLOT.begin(); it2!=NOPRISLOT.end(); ++it2)
	{
		fout<<(*it2)<<"  ";
	}
	fout<<"\n\n";
	for(IntSet::iterator it3=NOBACKPATH.begin(); it3!=NOBACKPATH.end(); ++it3)
	{
		fout<<(*it3)<<"  ";
	}
	fout<<"\n\n";
	for(IntSet::iterator it4=NOBACKSLOT.begin(); it4!=NOBACKSLOT.end(); ++it4)
	{
		fout<<(*it4)<<"  ";
	}
	fout<<"\n\n";
	for(IntSet::iterator it5=NOLAYER.begin(); it5!=NOLAYER.end(); ++it5)
	{
		fout<<(*it5)<<"  ";
	}
	fout<<"\n\n";
	for(IntSet::iterator it6=NOADDCON.begin(); it6!=NOADDCON.end(); ++it6)
	{
		fout<<(*it6)<<"  ";
	}
	fout<<"\n\n";
	*/
	fout << "\n\n----------------------END OF SIMULATION-----------------------\n";
}
