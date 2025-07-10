#include <iomanip>
#include <iostream>
#include <limits>
#include <algorithm>

#include "comm_brick.hpp"
#include "atom.hpp"
#include "read_vtp.hpp"
#include "integrate.hpp"
#include "read_dump.hpp"

using namespace PHAFD_NS;

/* dump file structure is going to be that each timestep is written to an individual (VTK)
   file, and then there is a .pvd file to display the timesteps collectively. Files are
   always going to be outputted in parallel. */
ReadVTP::ReadVTP(PHAFD *phafd) : Pointers(phafd) {

  Nfirst = -1;
  Nlast = -1;
  skip = -1;
  every = -1;
};



void ReadVTP::init(const std::vector<std::string> &v_line) {


  v_line_for_read_dump;


  filename = v_line.at(0);
  check_pvd_extension();
  size_t endslash = filename.find_last_of("\\/");

  if (endslash != std::string::npos) {
    directory = filename.substr(0,endslash+1);
  }

  auto it = v_line.begin()+1;
  while (it != v_line.end()) {
    if (*it == "first") {
      ++it;
      Nfirst = std::stoll(*it);
    } else if (*it == "last") {
      ++it;
      Nlast = std::stoll(*it);
    } else if (*it == "every") {
      ++it;
      every = std::stoll(*it);
    } else if (*it == "skip") {
      ++it;
      skip = std::stoll(*it);
    } else if (*it == "dump") {
      v_line_for_read_dump.assign(it,v_line.end());
      v_line_for_read_dump.at(0) = v_line_for_read_dump.at(1);
      v_line_for_read_dump.at(1) = "temporary";
      break;
    }
    ++it;
  }

  if (v_line_for_read_dump.size() == 0)
    throw std::runtime_error("ReadVTP requires dump keyword.");


  

  if (std::find(v_line.begin(),v_line.end(),"first") != v_line.end()
      && Nfirst == -1)
    throw std::runtime_error("ReadVTP requires dump keyword to be last.");
  
  if (std::find(v_line.begin(),v_line.end(),"last") != v_line.end()
      && Nlast == -1)
    throw std::runtime_error("ReadVTP requires dump keyword to be last.");

  if (std::find(v_line.begin(),v_line.end(),"every") != v_line.end()
      && every == -1)
    throw std::runtime_error("ReadVTP requires dump keyword to be last.");
  if (std::find(v_line.begin(),v_line.end(),"skip") != v_line.end()
      && skip == -1)
    throw std::runtime_error("ReadVTP requires dump keyword to be last.");  

  if (skip != -1 && every != -1)
    throw std::runtime_error("Cannot use both skip and every keyword in "
			     "ReadVTP.");
  
}

ReadVTP::~ReadVTP() {
}


void ReadVTP::read()
{

  read_preamble(filename,"<Collection>",
		vtpfile);


  int end_of_file = 0,stop;
  int global_err;

  int64_t skip_counter = 0;



  integrate->setup();
  while (true) {
    std::string pvtiname = read_vtp_for_pvti_name();

    
      
    
    if (pvtiname == "" && commbrick->me == 0)
      end_of_file = 1;

    MPI_Allreduce(&end_of_file,&stop,1,MPI_INT,MPI_SUM,world);

    if (stop)
      return;


    int64_t ptimestep = 0;

    if (commbrick->me == 0)
      ptimestep = check_pinstance_name(directory + pvtiname);

    MPI_Allreduce(&ptimestep,&global_err,1,MPI_INT64_T,MPI_MIN,world);

    if (global_err == -1)
      throw std::runtime_error(pvtiname
			       + std::string(" incorrect pvti naming")
			       + std::string(" format."));


    
  
    read_preamble(directory+pvtiname,"</PPointData>",
		  pvtifile);
  
    create_file_list();
    std::string vtiname = read_in_vti();

    int64_t timestep = check_instance_name(directory + vtiname);
    
    MPI_Bcast(&ptimestep,1,MPI_INT64_T,0,world);

    if (global_err == -1)
      throw std::runtime_error(pvtiname
			       + std::string(" incorrect vti naming")
			       + std::string(" format."));


    if (timestep != ptimestep)
      throw std::runtime_error("Different timesteps in pvti and vti"
			       " files.");


    std::vector<int64_t> cross_p_times(commbrick->nprocs);

    MPI_Allgather(&timestep,1,MPI_INT64_T,cross_p_times.data(),1,
		  MPI_INT64_T,world);


    for (int proc = 0; proc< commbrick->nprocs; proc++)
      if (timestep != cross_p_times.at(proc))
	throw std::runtime_error("Time steps don't match "
				 "across processors!");
    
    
    v_line_for_read_dump.at(1) = directory+vtiname;




    if (timestep < Nfirst) continue;
    else if (Nlast > 0 && timestep > Nlast) continue;
    else if (every > 0 && timestep % every != 0) continue;
    else if (skip > 0) {
      if (skip_counter % skip != 0) {
	skip_counter += 1;
	continue;
      }
      else skip_counter += 1;
    }


    
    if (commbrick->me == 0) {
      std::cout << "Rerunning " << pvtiname << " on timestep "
		<< timestep << std::endl;
    }

    ReadDump read_dump(phafd);
    read_dump.init(v_line_for_read_dump);
    read_dump.process_attributes();
    integrate->timestep = timestep;
    integrate->single_step(false);

    
  }

  integrate->finalise();
}

std::string::size_type ReadVTP::check_pvd_extension()
{
  std::string::size_type vstart;
  vstart = filename.find(".pvd");
  if (vstart == std::string::npos)
    throw std::runtime_error("invalid filename for ReadVTP");

  return vstart;
}


int64_t ReadVTP::check_instance_name(const std::string &vtiname)
{

  std::string::size_type vstart,vend;
  vstart = check_pvd_extension();

  std::string base_name = filename.substr(0,vstart);

  std::string expected;

  std::string shortened;

  expected = base_name + std::string("_p") + std::to_string(commbrick->me)
    + std::string("_");

  if (vtiname.substr(0,expected.length()) != expected) {
    std::cerr << vtiname.substr(0,expected.length())
	      << " but expected = "
	      << expected << std::endl;
    return -1;
  }

  shortened = vtiname.substr(expected.length());

  vstart = shortened.find(".vti");

  shortened = shortened.substr(0,vstart);


  int64_t step;

  try {
    step = std::stoll(shortened);
  } catch (std::invalid_argument &e) {
    std::cerr << shortened << " cannot be converted to int64_t"
	      << std::endl;

    return -1;
  }
  return step;
  
}

int64_t ReadVTP::check_pinstance_name(const std::string &pvtiname)
{

  std::string::size_type vstart,vend;
  vstart = check_pvd_extension();

  std::string base_name = filename.substr(0,vstart);

  std::string expected;

  std::string shortened;

  expected = base_name + std::string("_");

  if (pvtiname.substr(0,expected.length()) != expected) {
    std::cerr << pvtiname.substr(0,expected.length())
	      << " but expected = "
	      << expected << std::endl;
    return -1;
  }

  shortened = pvtiname.substr(expected.length());

  vstart = shortened.find(".pvti");

  shortened = shortened.substr(0,vstart);
  
  int64_t step;

  try {
    step = std::stoll(shortened);
  } catch (std::invalid_argument &e) {
    std::cerr << shortened << " cannot be converted to int64_t"
	      << std::endl;

    return -1;
  }

  return step;
  
}


  


      
  
  

void ReadVTP::read_preamble(const std::string &thename,
			    const std::string &keyword,
			    std::ifstream &file)
{

 std::string line;

 int localerr = 0,globalerr;
  if (commbrick->me == 0) {
    file = std::ifstream(thename);
    if (file.fail()) {
      std::cerr << "could not open file " << thename << std::endl;
      localerr = 1;
    }
  }
  MPI_Allreduce(&localerr,&globalerr,1,MPI_INT,MPI_SUM,world);
  
  if (globalerr)
    throw std::runtime_error("file not found.");
  
  
  if (commbrick->me == 0) {
    // push vtp file up to the line before the filenames are found
    std::string line;
    while (std::getline(file,line))
      if (line == keyword) break;
    
    if (file.eof())
      localerr = 1;
    
  }
  
  
  MPI_Allreduce(&localerr,&globalerr,1,MPI_INT,MPI_SUM,world);
  
  if (globalerr) 
    throw std::runtime_error("vtp file " + thename +
			     std::string(" is not in correct format."));
  
}


std::string ReadVTP::read_vtp_for_pvti_name()
{
  std::string line,subline = "";
  if (commbrick->me == 0) {
    std::getline(vtpfile,line);

    if (line == "</Collection>") return subline;
    
    std::string::size_type vstart,vend;

    vstart = line.find("file=\"");
    subline = line.substr(vstart+6);

    vend = subline.find("\"");

    subline = subline.substr(0,vend);

  }  

  
  return subline;
}

void ReadVTP::create_file_list()
{

  list_of_files.clear();
  int localerr=0,globalerr;

  if (commbrick->me == 0) {
    std::string line,subline;
    std::string::size_type vstart,vend;
    while (std::getline(pvtifile,line)) {

      if (line.find("<Piece") == std::string::npos)
	break;
      
      vstart = line.find("Source=\"");

      subline = line.substr(vstart+8);

      vend = subline.find("\"");
      subline = subline.substr(0,vend);
      list_of_files.push_back(subline);

    }
    if (pvtifile.eof())
      localerr = 1;
  }

  MPI_Allreduce(&localerr,&globalerr,1,MPI_INT,MPI_SUM,world);
  
  if (globalerr) 
    throw std::runtime_error("pvti file in " + filename +
			     std::string(" is not in correct format."));


    
  
  if (commbrick->me == 0) {
    if (list_of_files.size() != commbrick->nprocs)
      localerr = 1;
  }

  MPI_Allreduce(&localerr,&globalerr,1,MPI_INT,MPI_SUM,world);
  
  if (globalerr) 
    throw std::runtime_error("number of vti files specified in pvti "
			     "files must match number of processors.");

}


std::string ReadVTP::read_in_vti()
{


  std::string vtiname;
  if (commbrick->me == 0)  {
    int proc = 0;
    for (auto &st : list_of_files) {
      if (proc == 0)
	vtiname = st;
      else {
	MPI_Send(st.c_str(),st.length(),MPI_CHAR,proc,proc,world);
      }
      proc ++;
    }
  } else {
    MPI_Status status;
    MPI_Probe(0, commbrick->me, world, &status);
    int count;
    MPI_Get_count(&status,MPI_CHAR,&count);
    
    vtiname.resize(count);

    std::vector<char> buffer(vtiname.begin(), vtiname.end());
    MPI_Recv(&buffer[0],count,MPI_CHAR,0,commbrick->me,world,
	     &status);
    vtiname.assign(buffer.begin(), buffer.end());
    
  }


  return vtiname;
}
