#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <dirent.h>


int mark_done(char *tasks_done_dir, int myid)
{
	int rc = mkdir(tasks_done_dir, 0770);
	if (rc == -1 && errno != EEXIST) {
		fprintf(stderr, "failed to create done dir %s: %s\n", tasks_done_dir, strerror(errno));
		exit(1);
	}

	char task_done[1024];
	sprintf(task_done, "%s/task%d", tasks_done_dir, myid);
	int fds = open(task_done, O_CREAT | O_WRONLY, 0600);
	if (fds < 0)  {
		fprintf(stderr, "failed to create %s : %s\n", task_done, strerror(errno));
		exit(1);
	}
	close(fds);

	return 0;
}

int all_done(char *tasks_done_dir, int ntasks)
{
	DIR *dirp = opendir(tasks_done_dir);
	if (!dirp) {
		fprintf(stderr, "failed to open done dir %s: %s\n", tasks_done_dir, strerror(errno));
		exit(1);
	}

	struct dirent *de;

	int ntasks_done=0;
	errno = 0;
	for (de = readdir(dirp);
	     de;
	     de = readdir(dirp)) {
		ntasks_done++;
	}
	
	if (ntasks_done == ntasks+2) return 1; // remember . and ..

	if (errno) {
		fprintf(stderr, "failed to read done dir %s entry: %s\n", tasks_done_dir, strerror(errno));
		exit(1);
	}

	return 0;
}

struct task {
	int lon, lat;
};

struct task *read_tasklist(char *tasks_file)
{
	FILE *fd = fopen(tasks_file, "r");
	if (!fd)  {
		fprintf(stderr, "failed to open task file %s: %s\n", tasks_file, strerror(errno));
		exit(1);
	}

	int lon, lat, rc=0, count=0;
	while (rc != EOF) {	
		rc = fscanf(fd, "%d %d", &lon, &lat);
		printf("read %d %d\n", lon, lat);
		if (rc == 2) count++; 
	}
	printf("read %d tasks\n", count);
	
	struct task *tasklist = (task*) malloc((count+1)*sizeof(struct task));

	rewind(fd);

	for (int i=0; i<count;i++) {
		struct task *task = tasklist+i;
		rc = fscanf(fd, "%d %d", &task->lon, &task->lat);
	}

	fclose(fd);
	
	tasklist[count].lon = -500;

	return tasklist;
}


int main(int nargs, char **argv)
{
	// get my id
	int myid = atoi(argv[1]);
	int ntasks = atoi(argv[2]);
	char *tasks_file = argv[3];
	char tasklockfile[128];

	printf("args %d %d %s\n", myid, ntasks, tasks_file);
	struct task *task;
	struct task *screening_tasklist = read_tasklist(tasks_file);
	struct task *pairing_tasklist = screening_tasklist;
	// TODO read lists of tasks =>  screening_tasklist and pairing_tasklist

	char *bindir=getenv("PHES_BINDIR");
	if (!bindir) {
		fprintf(stderr, "PHES_BINDIR not set\n");
		exit(1);
	}

	char screening_cmd[1024];
	sprintf(screening_cmd, "%s/screening", bindir);
	char *screening_args = screening_cmd+strlen(screening_cmd);

	for (task = screening_tasklist; task->lon > -500; task++) {
		
		// form lockfile name
		sprintf(tasklockfile, "lockfiles/screening_task_%d_%d", task->lon, task->lat);

		printf("task %d trying %s\n", myid, tasklockfile);
		int fds = open(tasklockfile, O_CREAT | O_EXCL | O_WRONLY, 0600);
                if (fds < 0) {
                        if (errno == EEXIST) {
				// some other task has it
				continue;
                        } else {
				fprintf(stderr, "failed to to open lock file %s: %s\n", tasklockfile, strerror(errno));
				exit(1);
			}
		}

		
		printf("task %d got it\n", myid);
		// task is mine 
		sprintf(screening_args, " %d %d", task->lon, task->lat);
		int rc = system(screening_cmd);
		close(fds);
	}


	// all screening tasks have been allocated
	mark_done("screening_tasks_done", myid);

	while ( !all_done("screening_tasks_done", ntasks))
		sleep(5);

	printf("screening complete\n");

	char pairing_cmd[1024];
	sprintf(pairing_cmd, "%s/pairing", bindir);
	char *pairing_args = pairing_cmd+strlen(pairing_cmd);

	for (task = pairing_tasklist; task->lon > -500; task++) {
		
		// form lockfile name
		sprintf(tasklockfile, "lockfiles/pair_task_%d_%d", task->lon, task->lat);
		printf("task %d trying %s\n", myid, tasklockfile);

		int fds = open(tasklockfile, O_CREAT | O_EXCL | O_WRONLY, 0600);
                if (fds < 0) {
                        if (errno == EEXIST) {
				// some other task has it
				continue;
                        } else {
				fprintf(stderr, "failed to to open lock file %s: %s\n", tasklockfile, strerror(errno));
				exit(1);
			}
		}

		// task is mine 
		printf("task %d got it\n", myid);

		sprintf(pairing_args, " %d %d", task->lon, task->lat);
		int rc = system(pairing_cmd);
		close(fds);
	}

	
	// all pairing tasks have been allocated
	mark_done("pairing_tasks_done", myid);

	while ( !all_done("pairing_tasks_done", ntasks))
		sleep(5);

	printf("pairing complete\n");
}
