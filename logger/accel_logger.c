#include <stdio.h>
#include <inttypes.h>
#include <lcm/lcm.h>
#include "../lcmtypes/vicon_state_t.h"
#include "../lcmtypes/accel_t.h"
#include "../lcmtypes/command_t.h"

#include "../include/util.h"
#include <time.h>
#include <stdlib.h>

time_t rawtime;
struct tm* info;
char accel_fname[50];

void accel_handler(const lcm_recv_buf_t *rbuf, const char * channel, 
        const accel_t * msg, void * user) {

	FILE* logger;
	logger = fopen(accel_fname, "a");
	if(logger == NULL) {
		perror("File open error! \n");
		exit(0);
	}
	
	//timestamp, accel[0]~[2], angle[0]~[2], angular_vel[0]~[2]
	fprintf(logger, "%f, %f, %f\n", msg->accel[0], msg->accel[1], msg->accel[2]);
		fclose(logger);

}

int main(int argc, char* argv[]) {

	time(&rawtime);
	info = localtime(&rawtime);
	strftime(accel_fname, 50, "../data/accel_%d%H%M%S.txt", info);

	lcm_t* lcm = lcm_create(NULL);

	accel_t_subscribe(lcm, "acceleration_wrong", accel_handler, NULL);



	while(1) {
		lcm_handle(lcm);
	}

	lcm_destroy(lcm);
	return 0;
}
