#include <fcntl.h>
#include <dirent.h>

#include "phes_base.h"
#include "boost/filesystem.hpp"
namespace fs = boost::filesystem;

int display = true;

// Create a lockfile for this driver and process, returning the ID of the driver
int set_worker(string process){
	mkdir(convert_string(file_storage_location+"driver_files/"+process+"_workers"), 0770);
	int i = 0;
	while(true){
		string workerlockfile = file_storage_location+"driver_files/"+process+"_workers/"+to_string(i);
		int fds = open(convert_string(workerlockfile), O_CREAT | O_EXCL | O_WRONLY, 0600);
		
		if (!(fds < 0)) {
			mkdir(convert_string(file_storage_location+"driver_files/"+process+"_done_workers/"), 0770);
			mkdir(convert_string(file_storage_location+"driver_files/"+process+"_done_workers/"+to_string(i)), 0770);
			close(fds);
            return i;
		}
		close(fds);
		i++;
	}
}

// Deletes the lockfile with the given ID
void unset_worker(int id, string process){
	string workerlockfile = file_storage_location+"driver_files/"+process+"_done_workers/"+to_string(id)+"/done";
	int fds = open(convert_string(workerlockfile), O_CREAT | O_EXCL | O_WRONLY, 0600);
	close(fds);
	if (fds < 0)
        throw 1;
}

// Checks if all drivers have completed the process by checking if any files exist in the processes lockfile directory
bool all_done(string process)
{
	fs::path p(file_storage_location+"driver_files/"+process+"_done_workers/");
    std::vector<fs::directory_entry> v; // To save the file names in a vector.

    copy(fs::directory_iterator(p), fs::directory_iterator(), back_inserter(v));

    for ( std::vector<fs::directory_entry>::const_iterator it = v.begin(); it != v.end();  ++ it )
    {
        struct stat buffer;   
  		if(!(stat (((*it).path().string()+"/done").c_str(), &buffer) == 0))
  			return false; 
    }    

    return true;
}

// Reads a list of cells to process from the tasks_file (Eg. 148 -36)
vector<string> read_tasklist(char *tasks_file)
{
	ifstream fd(tasks_file);
	if (!fd)  {
		fprintf(stderr, "failed to open task file %s: %s\n", tasks_file, strerror(errno));
		exit(1);
	}
	vector<string> tasklist;
	string line;
	while(getline(fd, line)){
		line.erase(remove(line.begin(), line.end(), '\n'), line.end());
		line.erase(remove(line.begin(), line.end(), '\r'), line.end());
		tasklist.push_back(line);
	}
	printf("read %zu tasks\n", tasklist.size());
	fd.close();
	return tasklist;
}

// Reads a list of processes to complete from the processes_file (Eg. screening)
vector<string> read_processlist(char *processes_file)
{
	ifstream fd(processes_file);
	if (!fd)  {
		fprintf(stderr, "failed to open process file %s: %s\n", processes_file, strerror(errno));
		exit(1);
	}
	vector<string> processlist;
	string line;
	while(getline(fd, line)){
		line.erase(remove(line.begin(), line.end(), '\n'), line.end());
		line.erase(remove(line.begin(), line.end(), '\r'), line.end());
		processlist.push_back(line);
	}
	printf("read %zu processes\n", processlist.size());
	fd.close();
	return processlist;
}

void write_to_logfile(int id, string process, string message){
	ofstream logfile;
  	logfile.open(file_storage_location+"debug1/"+process+"_logfiles/"+process+"_"+to_string(id), ios_base::app);
  	logfile<<message;
  	logfile.close();
}

int main()
{
	parse_variables(convert_string("storage_location"));
	parse_variables(convert_string(file_storage_location+"variables"));

	vector<string> tasklist = read_tasklist(convert_string(file_storage_location+tasks_file));
	vector<string> processlist = read_processlist(convert_string(file_storage_location+processes_file));

	for (auto process : processlist) {
		mkdir(convert_string(file_storage_location+"debug1"), 0770);
		mkdir(convert_string(file_storage_location+"driver_files"), 0770);
		int id = set_worker(process);
		mkdir(convert_string(file_storage_location+"debug1/"+process+"_logfiles"), 0770);
		
		for(auto task : tasklist){
			mkdir(convert_string(file_storage_location+"driver_files/lockfiles"), 0770);
			string tasklockfile = file_storage_location+"driver_files/lockfiles/"+process+"_task_"+format_for_filename(task);
			int fds = open(convert_string(tasklockfile), O_CREAT | O_EXCL | O_WRONLY, 0600);
			if (fds < 0) {
                if (errno == EEXIST) {
					continue;
	            } else {
					fprintf(stderr, "failed to to open lock file %s: %s\n", convert_string(tasklockfile), strerror(errno));
					write_to_logfile(id, process, "Failed to open lock file "+tasklockfile+"\n");
					exit(1);
				}
			}
			close(fds);
			bool success = false;
			for(int i = 0; i<5; i++){
				if(!system(convert_string("./bin/"+process+" "+task))){
					success = true;
					break;
				}else{
					cout<<"Retrying command: ./bin/"+process+" "+task+"\n";
					write_to_logfile(id, process, "Retrying command: ./bin/"+process+" "+task+"\n");
					sleep(100);
				}
			}
			if(!success){
				cout<<"Problem running command: ./bin/"+process+" "+task+"\n";
				write_to_logfile(id, process, "Problem running command: ./bin/"+process+" "+task+"\n");
				// exit(1);
			}
			
		}
		unset_worker(id, process);
		write_to_logfile(id, process, "Done\n");
		while(!all_done(process)){
		 	sleep(10);
		}
	}
	printf("Done\n");
}
