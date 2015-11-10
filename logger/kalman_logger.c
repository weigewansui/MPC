#include <stdio.h>
#include <inttypes.h>
#include <lcm/lcm.h>
#include "../lcmtypes/vicon_state_t.h"
#include "../lcmtypes/kalman_v_t.h"

#include "../include/util.h"
#include <time.h>
#include <stdlib.h>

time_t rawtime;
struct tm* info;
char state_fname[50];
char kalman_fname[50];

void state_handler(const lcm_recv_buf_t *rbuf, const char * channel, 
        const vicon_state_t * msg, void * user) {

	FILE* logger;
	logger = fopen(state_fname, "a");
	if(logger == NULL) {
		perror("File open error! \n");
		exit(0);
	}
	
	//timestamp, accel[0]~[2], angle[0]~[2], angular_vel[0]~[2]
	fprintf(logger, "%" PRId64 ",%f,%f,%f,%f,%f,%f\n", msg->timestamp,
	 msg->position[0], msg->position[1], msg->position[2], msg->velocity[0], msg->velocity[1], msg->velocity[2]);
	fclose(logger);

}

void kalman_handler(const lcm_recv_buf_t *rbuf, const char * channel, 
        const kalman_v_t * msg, void * user) {

	FILE* logger;
	logger = fopen(kalman_fname, "a");
	if(logger == NULL) {
		perror("File open error! \n");
		exit(0);
	}
	
	//timestamp, accel[0]~[2], angle[0]~[2], angular_vel[0]~[2]
	fprintf(logger, "%" PRId64 ",%f,%f,%f,%f,%f,%f,%f,%f,%f\n",msg->timestamp, msg->pos[0], msg->pos[1], msg->pos[2], 
		msg->vel[0], msg->vel[1], msg->vel[2], msg->accel[0], msg->accel[1], msg->accel[2]);
	fclose(logger);

}

int main(int argc, char* argv[]) {

	time(&rawtime);
	info = localtime(&rawtime);
	strftime(state_fname, 50, "../data/vicon_state_true_%d%H%M%S.txt", info);
	strftime(kalman_fname, 50, "../data/kalman_state_%d%H%M%S.txt", info);

	printf("%s\n",state_fname);

	lcm_t* lcm = lcm_create(NULL);

	vicon_state_t_subscribe(lcm, "vicon_state_true", state_handler, NULL);
	kalman_v_t_subscribe(lcm, "kalman_v", kalman_handler, NULL);

	while(1) {
		lcm_handle(lcm);
	}

	lcm_destroy(lcm);
	return 0;
}