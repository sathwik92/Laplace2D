//=================================================================
//  File input_reader.cpp
//  Created:  March 17, 2017: Sathwik Bharadwaj
//  Modified: 
//  Last change: 
//=================================================================

// ****************************************************************
// This file reads in the input file and initializes the data 
//	structure. 
// ****************************************************************
// Functions within this file:
//	template <typename T> T readin(std::ifstream& file, const char* description, int ndebug);
//	int determine_current_line_number( std::ifstream& file);
//	void ignore_until ( std::ifstream& file, char c );
//=================================================================
//header file with all external libraries and user defined classes.
#include "global_params.h"
//=================================================================

/*!
 *
 * Returns the line in the file that the get pointer of the given 
 * std::ifstream is on.  \note The first line would be returned as one not zero.
 *
 * @param[in] file an std::ifstream specifying the file to be read from
 *
 *
 * This function is used by readin to determine the line number that input was read from.
 * 
 */
//===================================================================
int determine_current_line_number( std::ifstream& file)
{
	int pos = file.tellg();
	file.seekg(0, std::ios::beg);

	int linenum = 0;
	while(file.tellg() <= pos)
	{
		linenum++;
		ignore_until(file,'\n');
	}

	file.clear();
	file.seekg(pos, std::ios::beg);
	return linenum;
}

/*!
 *
 * Ignores all characters in the stream until the specified character
 * is reached.  The specified character is also ignored.
 *
 * @param[in] file an std::ifstream specifying the file to be read from
 * @param[in] c    a character specifying the character to stop at
 *
 *
 * This function is very similar to the ignore function that is a member of
 * the istream class, except that the istream function requires a maximum
 * search length to be specified, and once that many characters have been
 * searched through, the function will stop searching even if the specified
 * character has not been detected.  I could find no way to specify an 
 * unlimited search length, and so wrote this replacement.
 * 
 */

void ignore_until ( std::ifstream& file, char c )
{
	char ch;
	while ( file.get ( ch ) && ch != c ) ;
}
/*!
 *
 * @brief Allows data of a specified type to be readin from a file.  All 
 * whitespace is ignored except as separators between values, and a '#' 
 * character anywhere on a line will cause the rest of the line to be 
 * ignored. 
 *
 * @param[in]  file          an ifstream specifying the file to be read from
 * @param[out] data          a void* pointing to where the output is stored
 * @param[in]  description   a const char* containing the description of the value
 * @param[in]  type          a char containing the type of data to be read in.  For a list of recognized types see \ref RecognizedTypes.
 * @param[in]  units         a const char* specifying what the units (if any) of the input are.
 * @param[in]  ndebug        a int specifying the debug output level.  For a list of recognized debug levels see \ref DebugLevels. 
 *
 * \section WhitespaceAndCommentHandling Whitespace and Comment Handling
 * Whitespace is completely ignored except as a separator for data values.  
 * Because of this, there is no way to load a string with whitespace in it,
 * or a char containing whitespace.  Also the " and ' characters are treated 
 * as regular characters.
 *
 * If at any Point in the file a '#' character is encountered, the rest of
 * the line will be ignored.  As a result, there is no way to load the '#'
 * character in a string or in a character.
 *
 * 
 * \section Warnings Possible Error and Warning Messages
 * \par
 * If the end of the file is reached while reading for input, an exception 
 * is thrown and a message will be printed specifying the description and 
 * type of the input that caused the error.  
 *
 * \par
 * If the specified type was not recognized an exception is thrown and a 
 * message is printed specifying the description and type of the input, and 
 * the list of recognized inputs.For a list of recognized types data types see 
 * \ref RecognizedTypes.
 *
 * \par
 * For warnings regarding conversion, see \ref DataConversion.   
 * 
 * \section DataConversion Description of Data Conversion
 * All conversions are done using stringstreams.  Data is read from the file as
 * a string.  This string is pushed into a stringstream, and then data from the
 * stringstream is extracted into the output with the specified type (for a list
 * of recognized types see \ref RecognizedTypes).  If after doing this the
 * stringstream's failbit or badbit is set, then an error is raised, and a
 * message printed out stating that conversion was not possible and specifying
 * the line where the input was read from, the type and description of  the
 * expected input, and the input string that was readin.  If after doing this the
 * stringstream's failbit and badbit are not set, but neither is the eofbit, then
 * a warning is printed stating that some of the data that was readin could not
 * be converted and specifying the actual read in string and the type that it
 * tried to convert it to.
 */
	template <typename T> 
T readin(std::ifstream& file, const char* description, int ndebug)
{
	std::stringstream s;
	std::string str;
	// Read next input, skipping white space
	file >> std::skipws >> str;
	T ansi;
	int i = 0;

	// If this input is a comment
	while(str[0] == '#')
	{
		// Ignore the rest of the line
		ignore_until ( file, '\n');
		// and continue reading
		file >> std::skipws >> str;
		if(file.eof())
		{
			std::cout << "ERROR reached end of file while scanning for input" << description << "." << std::endl;
			throw(1);
		}
	}
	if(file.eof())
	{
		std::cout << "ERROR reached end of file while scanning for input" << description << "." << std::endl;
		throw(1);
	}
	// Search for any #'s that may have been read in accidentaly
	// for example if a file had
	//
	// 5#This is a comment
	//
	// then "5#This" would be be read in, and we only want the 5 to
	// be read in.  So, we put the "#This" back into the stream and 
	// remove it from str
	//
	// First, we determine where the first '#' is, if there is one
	int comment_location = str.length();
	//determine the location
	i = 0;
	for(auto svalue: str)
	{
		if(svalue == '#')
		{
			comment_location = i;
			break;
		}
		i++;
	}
	// If comment_location==str.length(), then that means this string
	if(comment_location != str.length())
	{
		// Loop over the string backwards until we reach
		// the location of the first '#' and put back each
		// character into the stream
		for(i=str.length()-1; i >= comment_location; i--)
		{
			file.putback(str[i]);
		}
		// Then remove the comments from the string
		str.erase(str.begin()+comment_location,str.end());
	}

	// Put the string into the steam to prepare for conversion
	s << str;
	// Attempt a conversion
	s >> ansi;
	// Detect and Handle Errors and Warnings
	if(s.fail()) 
	{
		std::cout << "ERROR on line ";
		std::cout.fill('0');
		std::cout.width(4);
		std::cout << determine_current_line_number(file)
			<< " expected integer for input " << description << 
			", but input was \"" << str << "\". Conversion "
			"was not possible" << std::endl;
		throw(1);
	}
	if(ndebug > 0)
	{
		std::cout << "Line:\t" << determine_current_line_number(file) <<"\t"<< description << " " << ansi << "\n";
	}

	if(!s.eof() and ndebug > 0) 
	{
		std::cout << "----> expected integer input, but input was \"" << str <<  
			"\". However, a partial conversion was possible." << std::endl;
	}
	else if(ndebug > 0)
		if(!s.eof() and ndebug > 1) 
			throw(1);

	return ansi;
}

/*
 * This function is used to read the parameters from the input file in.  
 */
void input_reader(std::ifstream &in, data &dat)
{

	//=================================================================
	int i; //dummy variable
	//=================================================================

	dat.ndebug = readin<int>(in,"Debug Level:", 1);       // 'I' - integer

	dat.xmin = readin<double>(in, "Minimum X coordinate:", dat.ndebug); 

	dat.xmax = readin<double>(in, "Maximum X coordinate:", dat.ndebug);

	dat.ymin = readin<double>(in, "Minimum Y coordinate:" , dat.ndebug);

	dat.ymax = readin<double>(in, "Maximum Y coordinate", dat.ndebug);

	dat.ndzx = readin<int>(in,"Data points along x axis for plotting:", dat.ndebug);

	dat.ndzy = readin<int>(in,"Data points along y axis for plotting:", dat.ndebug);
	//=================================================================

	// dat.dzx, dat.dzy determines the step between two interpolation points.

	dat.dzx = (dat.xmax-dat.xmin)/(double)dat.ndzx;
	dat.dzy = (dat.ymax-dat.ymin)/(double)dat.ndzy;

	//=================================================================

	dat.node_elem = readin<int>(in, "Nodes per element:", dat.ndebug);

	dat.ndof = readin<int>(in,"Degrees of Freedom:", dat.ndebug);

	dat.nband = readin<int>(in, "Number of Bands:", dat.ndebug);

	dat.ngaus = readin<int>(in,"Order  of Gauss quadrature:", dat.ndebug);

	//=================================================================

	// Solver information:

	//=================================================================

	dat.rtol = readin<double>(in, "Relative convergence tolerance:", dat.ndebug);

	dat.maxit = readin<int>(in, "Maximum allowed iterations:", dat.ndebug);

	dat.atol = readin<double>(in, "Absolute convergence tolerance:", dat.ndebug);

	dat.divtol = readin<double>(in, "Divergence tolerance:", dat.ndebug);

	dat.solvertype = readin<std::string>(in, "Solver type:", dat.ndebug);

	dat.output_files_path = readin<std::string>(in, "Output_files_path:", dat.ndebug);   

	//=================================================================
	// Now set up the FEM parameters:
	// The max number of Gauss points in the file utiltrig.cpp
	//=================================================================
	dat.xigaus = new double[dat.ngaus];	
	dat.etagaus = new double[dat.ngaus];
	dat.wgaus = new double[dat.ngaus];

	trigauss(dat.ngaus, dat.xigaus, dat.etagaus, dat.wgaus);

	//=========================================

	dat.nglobal  =  dat.ngnodes*dat.ndof*dat.nband;    // Number of rows in global matrix 
	dat.nele_mat =  dat.node_elem*dat.ndof*dat.nband;  // Number of rows in element matrix
	dat.ndeg     =  dat.node_elem*dat.ndof;      // Degree of interpolation polynomials - not really degree of polynomial
						     // more like # of interpolation
						     // polynomials

	//=================================================================
	// open file to output the actual data used: this is for verification
	//=================================================================
	char outputpath[50];
	strncpy(outputpath, dat.output_files_path.c_str(), sizeof(outputpath));
	outputpath[sizeof(outputpath) - 1] = 0;
	char buffer[200];
	sprintf(buffer,"%sinput_params.out", outputpath);
	std::ofstream fout(buffer);
	fout <<"====================================" << std::endl;
	fout <<"============================" << std::endl;
	fout <<"Example FEM calculation: Laplace equations in 2D" << std::endl;
	fout <<"      Sathwik Bharadwaj     " << std::endl;
	fout <<"============================" << std::endl;
	fout <<"Size of the global matrix and vector: "<< dat.nglobal <<std::endl;
	fout <<"Global number of nodes     : " << dat.ngnodes << std::endl<< std::endl;
	fout <<"Number of elements         : " << dat.nelem <<std::endl;

	fout <<"Number of nodes per element: " << dat.node_elem << std::endl;
	fout <<"Degrees of Freedom per node: " << dat.ndof << std::endl;
	fout <<"Size of element matrix   : " << dat.nele_mat << std::endl;  
	fout <<"============================" << std::endl;
	fout <<"====================================" << std::endl;
	fout.close();
	return;
} // End of subroutine input_reader


// ==============================
// Mesh reader 
// ==============================
PetscErrorCode mesh_input(data &dat) 
{
	PetscErrorCode ierr;//Petsc Error Code
	std::ifstream ine,inn,inbn,inbe; //input files
	ine.open("../mesh/elem.dat"); //Element data. Found in mesh/ file
	inn.open("../mesh/node.dat"); //Node data. Found in mesh/ file
	inbn.open("../mesh/bnode.dat"); //Boundary Node data. Found in mesh/ file
	inbe.open("../mesh/belem.dat"); //Boundary Element data. Found in mesh/ file
					// opening mesh output files
	int i,j; // dummy variables
	int itrash, bnumber;
	if(!ine)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD, "ERROR: The file elem.dat was not found");CHKERRQ(ierr);
		exit(1);
	}

	dat.nelem = readin<int>(ine, "Number of elements:",dat.ndebug); //obtaining #of elements used in readin loop later
	dat.elem = new int*[dat.nelem];  //declaring dynamic arrays - 1st part
	dat.material = new int[dat.nelem];

	for(i=0; i< dat.nelem; i++)
	{
		dat.elem[i] = new int[3];  //declaring dynamic arrays - 2nd part
		j = readin<int>(ine, "",0);
		dat.elem[i][0] = readin<int>(ine, "",0);  //the three different nodes of the triangular elements
		dat.elem[i][1] = readin<int>(ine, "",0);  //indexing is made in a counterclockwise way
		dat.elem[i][2] = readin<int>(ine, "",0);
		dat.material[i] = readin<int>(ine,"",0);
		itrash = readin<int>(ine, "", 0); //neighbouring nodes from the elem.mesh are not required. 
		itrash = readin<int>(ine, "", 0);
		itrash = readin<int>(ine, "", 0);
	}
	ine.close();  //closing std::ifstream ine
	if(!inn)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD,"ERROR: The file node.dat was not found.\n");CHKERRQ(ierr);
		exit(1);
	}
	dat.ngnodes = readin<int>(inn,"Number of Nodes:",dat.ndebug);
	dat.node = new double*[dat.ngnodes];

	for(i=0; i< dat.ngnodes; i++)
	{
		dat.node[i] = new double[2];
		j = readin<int>(inn, "",0);
		dat.node[i][0] = readin<double>(inn, "",0);  //x-coordinate of the node
		dat.node[i][1] = readin<double>(inn, "",0);  //y-coordinate of the node
	}
	inn.close(); //closing std::ifstream inn
	if(!inbn)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD,  "ERROR: The file bnode.dat was not found.\n");CHKERRQ(ierr);
		exit(1);
	}

	dat.nbnode = readin<int>(inbn,"Number of Boundary Nodes:",dat.ndebug);
	dat.nboundary = readin<int>(inbn,"Number of Boundaries:", dat.ndebug);
	dat.bnode = new double*[dat.nbnode];
	dat.bvalue = new double*[dat.nbnode];
	dat.btob = new int[dat.nbnode];
	for(i=0; i<dat.nbnode; i++)
	{
		dat.bnode[i] = new double[2];
		dat.bvalue[i] = new double[6];
		//     dat.bnum[i]  = new double[1];
		j = readin<int>(inbn,"",0); 
		dat.bnode[i][0] = readin<double>(inbn, "", 0);  //corresponds to x-coordinate of the bnode
		dat.bnode[i][1] = readin<double>(inbn, "", 0); //corresponds to the y-coordinate of the bnode
		bnumber = readin<int>(inbn, "" , 0);//bnumber will provide the boundary number 
						    // boundary conditions
						    //Lower horizontal line
		if(bnumber ==1)
		{
			dat.bvalue[i][0] = 0; //bvalue[i][0] is the "function" value
			dat.bvalue[i][1] = 0; //bvalue[i][1] - derivative with respect to x
			dat.bvalue[i][2] = 99; //bvalue[i][2] - derivative with respect to y. Normal derivatives are left floating.    
			dat.bvalue[i][3] = 0; //bvalue[i][3] - 2nd derivative with respect to x,x        
			dat.bvalue[i][4] = 99; //bvalue[i][4] - 2nd derivative with respect to x,y        
			dat.bvalue[i][5] = 99; //bvalue[i][5] - 2nd derivative with respect to y,y    
			dat.btob[i] = 0; //type of boundary condtion. Here its dirchlet.
		}
		// right vertical line
		else if(bnumber ==2)
		{
			dat.bvalue[i][0] = 0; //bvalue[i][0] is the "function" value
			dat.bvalue[i][1] = 99; //bvalue[i][1] - derivative with respect to x
			dat.bvalue[i][2] = 0; //bvalue[i][2] - derivative with respect to y. Normal derivatives are left floating.    
			dat.bvalue[i][3] = 99; //bvalue[i][3] - 2nd derivative with respect to x,x        
			dat.bvalue[i][4] = 99; //bvalue[i][4] - 2nd derivative with respect to x,y        
			dat.bvalue[i][5] = 0; //bvalue[i][5] - 2nd derivative with respect to y,y    
			dat.btob[i] = 0; //type of boundary condtion. Here its dirchlet.
		}
		//upper horizontal line
		else if(bnumber ==3)
		{
			dat.bvalue[i][0] = 0; //bvalue[i][0] is the "function" value
			dat.bvalue[i][1] = 0; //bvalue[i][1] - derivative with respect to x
			dat.bvalue[i][2] = 99; //bvalue[i][2] - derivative with respect to y. Normal derivatives are left floating.    
			dat.bvalue[i][3] = 0; //bvalue[i][3] - 2nd derivative with respect to x,x        
			dat.bvalue[i][4] = 99; //bvalue[i][4] - 2nd derivative with respect to x,y        
			dat.bvalue[i][5] = 99; //bvalue[i][5] - 2nd derivative with respect to y,y    
			dat.btob[i] = 0; //type of boundary condtion. Here its dirchlet.
		}
		// right vertical line
		else// if(bnumber ==4)
		{
			dat.bvalue[i][0] =10*sin(PI*(1.0-((dat.bnode[i][1])/20.0))); //bvalue[i][0] is the "function" value
			if(i ==0)
				dat.bvalue[i][1] = 0; //x derivative is zero if i=0 to be consistent with previous line
			else
				dat.bvalue[i][1] = 99;
			//////////////////////////////////////
			dat.bvalue[i][2] = 10*(PI/20.0)*cos(PI*(1.0-(dat.bnode[i][1]/20.0))); //bvalue[i][2] - derivative with respect to y. Normal derivatives are left floating.    
			if (i ==0)
				dat.bvalue[i][3] = 0; //bvalue[i][3] - 2nd derivative with respect to x,x
			else
				dat.bvalue[i][3] = 99;// to be consistent with the previous node
			dat.bvalue[i][4] = 99; //bvalue[i][4] - 2nd derivative with respect to x,y        
			dat.bvalue[i][5] = -10*(PI/20.0)*(PI/20.0)*sin(PI*(1-(dat.bnode[i][1]/20.0))); //bvalue[i][5] - 2nd derivative with respect to y,y    
			dat.btob[i] = 0; //type of boundary condtion. Here its dirchlet.
		}
	}
	inbn.close(); //closing std::ifstream inbn

	if(!inbe)
	{
		ierr = PetscPrintf(PETSC_COMM_WORLD,"ERROR: The file belem.dat was not found.\n");CHKERRQ(ierr);
		exit(1);
	}
	dat.nbelem = readin<int>(inbe, "No. of Boundary Elements",  dat.ndebug);
	dat.belem = new int*[dat.nbelem];
	dat.bmaterial = new int*[dat.nbelem];
	for(i=0; i<dat.nbelem; i++){
		dat.belem[i] = new int[2];
		dat.bmaterial[i] = new int[2];
		j = readin<int>(inbe, "",0);
		dat.belem[i][0] = readin<int>(inbe, "",0);  //first node from which the element starts
		dat.belem[i][1] = readin<int>(inbe, "", 0);  //second node at which the element ends
		dat.bmaterial[i][0] = readin<int>(inbe, "",0); //element to the right of the element(counterclockwise)
		dat.bmaterial[i][1] = readin<int>(inbe, "",0); //element to the left - not a necessary thing. We can even put into trash.
	}
	inbe.close(); //closing std::ifstream inbe

	// ==========================================================================
	// the next section test the mesh_input function by recreating the mesh 
	// inputs from the dat struct
	// ==========================================================================

	std::ofstream felem;
	felem.open("../mesh/elem2.dat");
	for(i=0;i<dat.nelem;i++){
		felem<<i<<"\t"<<dat.elem[i][0]<<"\t"<<dat.elem[i][1]<<"\t"<<dat.elem[i][2]<<
			"\t"<<dat.material[i]<<std::endl;
	}
	felem.close();

	std::ofstream fnode;
	fnode.open("../mesh/node2.dat");

	for(i=0;i<dat.ngnodes;i++){
		fnode<<i<<"\t"<<dat.node[i][0]<<"\t"<<dat.node[i][1]<<std::endl;
	}
	fnode.close(); 


	std::ofstream fbnode;
	fbnode.open("../mesh/bnode2.dat");

	for(i=0;i<dat.nbnode;i++){
		fbnode<<i<<"\t"<<dat.bnode[i][0]<<"\t"<<dat.bnode[i][1]<<"\t"<<
			dat.bvalue[i][0]<<"\t"<<dat.bvalue[i][1]<<"\t"<<dat.bvalue[i][2]<<"\t"<<
			dat.bvalue[i][3]<<"\t"<<dat.bvalue[i][4]<<"\t"<<dat.bvalue[i][5]<<"\t"<<
			dat.btob[i]<<std::endl;
	}
	fbnode.close(); 

	std::ofstream fbelem;
	fbelem.open("../mesh/belem2.dat");
	for(i=0;i<dat.nbelem;i++){
		fbelem<<i<<"\t"<<dat.belem[i][0]<<"\t"<<dat.belem[i][1]
			<<"\t"<<dat.bmaterial[i][0]<<"\t"<<dat.bmaterial[i][1]<<std::endl;
	}
	fbelem.close(); 
	return ierr;
}
//===================================================================
