#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <dirent.h>

//Being risky (for devel only)
#include <bits/stdc++.h>
using namespace std;

#include "phes_base.h"

int display = false;

int set_worker(string process){
	mkdir(convert_string("driver_files/"+process+"_workers"), 0770);
	int i = 0;
	while(true){
		string workerlockfile = "driver_files/"+process+"_workers/"+to_string(i);
		int fds = open(convert_string(workerlockfile), O_CREAT | O_EXCL | O_WRONLY, 0600);
		if (!(fds < 0)) {
            return i;
		}
		i++;
	}
}

void unset_worker(int id, string process){
	string workerlockfile = "driver_files/"+process+"_workers/"+to_string(id);
	if(remove(convert_string(workerlockfile))!=0){
		printf("Problem removing worker: %d\n", id);
	}
}

bool all_done(string process)
{
	DIR *dir = opendir(convert_string("driver_files/"+process+"_workers"));
	if (!dir) {
		fprintf(stderr, "failed to open done dir %s: %s\n", convert_string("driver_files/"+process+"_workers"), strerror(errno));
		exit(1);
	}

	int n=0;
	while(readdir(dir)!=NULL) n++;
	if(n<=2){
		return true;
	}
	return false;
}


vector<GridSquare> read_tasklist(char *tasks_file)
{
	FILE *fd = fopen(tasks_file, "r");
	if (!fd)  {
		fprintf(stderr, "failed to open task file %s: %s\n", tasks_file, strerror(errno));
		exit(1);
	}

	vector<GridSquare> tasklist;

	int lon, lat, rc=0;
	while (rc != EOF) {	
		rc = fscanf(fd, "%d %d", &lon, &lat);
		tasklist.push_back(GridSquare_init(lat, lon));
	}
	printf("read %zu tasks\n", tasklist.size());
	fclose(fd);
	return tasklist;
}

vector<string> read_processlist(char *processes_file)
{
	ifstream fd(processes_file);
	if (!fd)  {
		fprintf(stderr, "failed to open process file %s: %s\n", processes_file, strerror(errno));
		exit(1);
	}

	vector<string> processlist;

	for(string line; getline(fd, line);){
		processlist.push_back(line);
	}
	printf("read %zu processes\n", processlist.size());
	fd.close();
	return processlist;
}


int main(int nargs, char **argv)
{
	char *tasks_file = argv[1];
	char *process_file = argv[2];

	printf("args %s %s\n", tasks_file, process_file);

	vector<GridSquare> tasklist = read_tasklist(tasks_file);
	vector<string> processlist = read_processlist(process_file);

	char *temp=getenv("PHES_BINDIR");
	if (!temp) {
		fprintf(stderr, "PHES_BINDIR not set\n");
		exit(1);
	}
	string bindir(temp);

	for (auto process : processlist) {
		int id = set_worker(process);
		for(auto task : tasklist){
			mkdir(convert_string("driver_files/lockfiles"), 0770);
			string tasklockfile = "driver_files/lockfiles/"+process+"_task_"+to_string(task.lon)+"_"+to_string(task.lat);
			int fds = open(convert_string(tasklockfile), O_CREAT | O_EXCL | O_WRONLY, 0600);
			if (fds < 0) {
                if (errno == EEXIST) {
					// some other program has it
					continue;
	            } else {
					fprintf(stderr, "failed to to open lock file %s: %s\n", convert_string(tasklockfile), strerror(errno));
					exit(1);
				}
			}
			if(system(convert_string(bindir+"/"+process+" "+to_string(task.lon)+" "+to_string(task.lat)))){
				printf("Problem running command\n");
				exit(1);
			}
		}
		unset_worker(id, process);
		while(!all_done(process)){
		 	sleep(5);
		}
	}
	printf("\nDone\n");
}
