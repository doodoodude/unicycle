cmake_minimum_required(VERSION 2.8.3)
project(unicycle)

find_package(catkin REQUIRED COMPONENTS
  controller_manager
  gazebo_ros
  gazebo_ros_control
  joint_state_controller
  position_controllers
  robot_state_publisher
  roscpp
  rviz
  xacro
)

catkin_package()

install(DIRECTORY launch
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})

install(DIRECTORY worlds
  DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION})

include_directories(include
  ${catkin_INCLUDE_DIRS}
)

add_executable(unicycle_controller src/unicycle_controller.cpp)
add_executable(stand_controller src/stand_controller.cpp)

target_link_libraries(unicycle_controller ${catkin_LIBRARIES})
target_link_libraries(stand_controller ${catkin_LIBRARIES})


