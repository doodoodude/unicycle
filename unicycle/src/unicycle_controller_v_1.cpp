#include "ros/ros.h"
#include "math.h"
#include "std_msgs/Float64.h"
#include "sensor_msgs/Imu.h"
#include "sensor_msgs/JointState.h"
#include <sstream>
#include <tf/transform_datatypes.h>

class Unicycle_controller
{
private:
  ros::NodeHandle n_;   

  ros::Subscriber imu_data_sub_;
  ros::Subscriber joint_states_sub_;

  ros::Publisher wheel_; 
  ros::Publisher disc_; 

  ros::Publisher roll_pub_; 

public: 

  double imu_ang_vel_x, imu_ang_vel_y, imu_ang_vel_z,
         imu_ang_acc_x, imu_ang_acc_y, imu_ang_acc_z,
         imu_lin_acc_x, imu_lin_acc_y, imu_lin_acc_z,
         imu_roll, imu_pitch, imu_yaw,

         err_integral,

         wheel_ang, disc_ang,
         wheel_vel, disc_vel,

         roll_min, roll_max, roll_abs,
         pitch_max,pitch_abs;            

  bool STOP;

  std_msgs::Float64 wheel_com, disc_com, 
                    roll_pub_data;

  Unicycle_controller(){ 
    imu_data_sub_ = n_.subscribe("/unicycle/imu_data",10,&Unicycle_controller::ImuDataReceive,this);
    joint_states_sub_ = n_.subscribe("/unicycle/joint_states",10,&Unicycle_controller::JointStatesObserve,this);

    wheel_ = n_.advertise<std_msgs::Float64>("/unicycle/wheel_controller/command",1);
    disc_ = n_.advertise<std_msgs::Float64>("/unicycle/disc_controller/command",1);

    roll_pub_ = n_.advertise<std_msgs::Float64>("/unicycle/roll",1); 

    roll_max = 1.4;
    roll_min = 0;
    pitch_max = 1.57;

  }

  void JointStatesObserve(const sensor_msgs::JointState::ConstPtr& joint_state){

    disc_ang = joint_state->position[0];
    wheel_ang  = joint_state->position[1];

    disc_vel = joint_state->velocity[0];
    wheel_vel  = joint_state->velocity[1];

  }

  void ImuDataReceive(const sensor_msgs::Imu::ConstPtr& imu_data){

    //imu_lin_acc_x = imu_data->linear_acceleration.x;
    //imu_lin_acc_y = imu_data->linear_acceleration.y;

    //imu_ang_acc_x = (imu_data->angular_velocity.x - imu_ang_vel_x)*100; //cм.внизу
    //imu_ang_acc_y = (imu_data->angular_velocity.y - imu_ang_vel_y)*100;

    imu_ang_vel_x = imu_data->angular_velocity.x;   
    imu_ang_vel_y = imu_data->angular_velocity.y;

    tf::Quaternion q(imu_data->orientation.x, imu_data->orientation.y, imu_data->orientation.z, imu_data->orientation.w);
    tf::Matrix3x3 orientation_matrix(q);
    orientation_matrix.getRPY(imu_roll, imu_pitch, imu_yaw);

    roll_abs = fabs(imu_roll);
    pitch_abs = fabs(imu_pitch);

//проверяем, не упал ли робот
    if(pitch_abs<pitch_max&&roll_abs<roll_max){
      STOP = false;
    }else{
      STOP = true;
    }

    roll_pub_data.data = imu_roll;
    roll_pub_.publish(roll_pub_data);

    pitch_stabilisation();
    roll_stabilisation();

  }

  void pitch_stabilisation(){

    err_integral = err_integral + imu_pitch/200;

    if(!STOP){
//если робот не упал

      wheel_com.data = (-15)*imu_pitch + (-1)*imu_ang_vel_y + (0)*err_integral;

    }else{

//иначе тормозим колёса
      wheel_com.data = (-0.05)*wheel_vel;
    }

    wheel_.publish(wheel_com);

  }

  void roll_stabilisation(){


    if(!STOP){
//если робот не упал

      if(roll_abs>roll_min){
  //если робот НЕ около положения равновесия

        disc_com.data = ((0)*disc_ang + (0)*disc_vel + (14.1908)*imu_roll + (2.3233)*imu_ang_vel_x)*2;

      }else{
  //иначе ничо не делаем

        disc_com.data = 0;
        
      }

    }else{
//если робот не упал
        disc_com.data = (-0.01)*disc_vel;
    }

    disc_.publish(disc_com);

  }

};


int main(int argc, char **argv)
{
  ros::init(argc, argv, "unicycle_controller");

  Unicycle_controller unicycle_controller;   

  ROS_INFO("Unicycle controller is running!");

  ros::spin();

/*
  ros::Rate loop_rate(200);                
  int count = 0;
  while (ros::ok())
    {
      ros::spin();
      loop_rate.sleep();
      ++count;
    }
*/
  
  return 0;
 }

