#include "ros/ros.h"
#include "math.h"
#include "std_msgs/Float64.h"
#include "sensor_msgs/JointState.h"
#include <sstream>

class stand_controller
{
public:

  ros::NodeHandle n;   

  double s, 
  post_p_ang, post_r_ang,
  stand_x_pos, stand_y_pos,
  post_p_amp, post_r_amp,
  stand_x_amp, stand_y_amp,
  delay,
  stop_delay;

  ros::Subscriber joint_state_sub;

  ros::Publisher post_p_pub,post_r_pub,stand_x_pub,stand_y_pub; 
  
  std_msgs::Float64 pub_data; 

  stand_controller(){

    post_p_amp = 0.1; post_r_amp = 0.1; stand_x_amp = 1.4; stand_y_amp = 0.4;

    delay = 3;
    stop_delay = 14;

    joint_state_sub = n.subscribe("/stand/joint_states",200,&stand_controller::joint_states_receive,this);

    post_p_pub = n.advertise<std_msgs::Float64>("/stand/post_p_controller/command",1);
    post_r_pub = n.advertise<std_msgs::Float64>("/stand/post_r_controller/command",1);
    stand_x_pub = n.advertise<std_msgs::Float64>("/stand/stand_x_controller/command",1);
    stand_y_pub = n.advertise<std_msgs::Float64>("/stand/stand_y_controller/command",1);

  }

  void joint_states_receive(const sensor_msgs::JointState::ConstPtr& joint_state){

    post_p_ang = joint_state->position[0];
    post_r_ang = joint_state->position[1];
    stand_x_pos = joint_state->position[2];
    stand_y_pos = joint_state->position[3];

    do_something();

  }

  void do_something(){

    s = ros::Time::now().toSec();

    if(s>=delay&&s<stop_delay){

      s = s - delay;

      pub_data.data = post_p_amp*sin(s*0.2);
      post_p_pub.publish(pub_data);

      pub_data.data = -post_r_amp*sin(s*0.3);
      post_r_pub.publish(pub_data);
      
      pub_data.data = -stand_x_amp*cos(s*0.2)+stand_x_amp;
      stand_x_pub.publish(pub_data);

      pub_data.data = stand_y_amp*sin(s*0.3);
      stand_y_pub.publish(pub_data);
      
    }

  }

}; 


int main(int argc, char **argv)
{
  ros::init(argc, argv, "stand_controller");

  stand_controller stand_cntrl;

  ROS_INFO("Stand controller is running!");

  ros::spin();
  
  return 0;
 }

