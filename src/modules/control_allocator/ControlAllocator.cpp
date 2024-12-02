/**************************
 *
 *   Copyright (c) 2013-2019 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 **************************/

/**
 * @file ControlAllocator.cpp
 *
 * Control allocator.
 *
 * @author Julien Lecoeur <julien.lecoeur@gmail.com>
 */
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <cmath>
#include "ControlAllocator.hpp"

#include <drivers/drv_hrt.h>
#include <circuit_breaker/circuit_breaker.h>
#include <mathlib/math/Limits.hpp>
#include <mathlib/math/Functions.hpp>


//ideaforge_team86
//declaring things failure_flag.h in the code
#include <uORB/topics/failure_flag.h> 
//control_allocator // printing
#include <px4_platform_common/events.h>
#include <px4_platform_common/log.h>
#include <matrix/math.hpp>
#include <uORB/topics/vehicle_attitude.h>
#include <uORB/topics/vehicle_angular_velocity.h>
#include <uORB/topics/vehicle_local_position.h>
#include <uORB/topics/sensor_gps.h>
#include <uORB/topics/vehicle_status.h>
#include <uORB/Subscription.hpp>
#include <px4_platform_common/defines.h>
#include <uORB/topics/vehicle_attitude_setpoint.h>

using namespace matrix;
using namespace time_literals;
//declaration of variables

//for failure detection
bool fail_change;
int faillM_indx=3;
bool check=false;
bool geo;

//for control algo
float roll, pitch, yaw;
float roll_rate, pitch_rate, yaw_rate;
float altitude,position_x,position_y;
float velocity_x,velocity_y,velocity_z;
float desired_roll, desired_pitch, desired_yaw, thrust_setpoint;
double latitude, longitude;
Vector3f torque_output;
float thrust_output;
Vector4f motor_outputs;
float pi = 3.14;
//control parammeters for geometric controller-1
static float roll_sp, pitch_sp, yaw_sp;
static bool calibrate;
constexpr float KP_ATTITUDE = 5.0f;    // Proportional gain for attitude
constexpr float KD_ATTITUDE = 0.5f;    // Derivative gain for angular rates
constexpr float MAX_THRUST = 1.0f;     // Max thrust for each motor
constexpr float MIN_THRUST = 0.0f; 
constexpr float L = 0.2f;              // Lever arm length (meters)
constexpr float D = 0.1f; 
/////////////////////////////************************************************************//////////////////////////////////////////// */
//main control parameters for geometric control
const double mass = 0.61;      // Mass of the drone (kg)
const double k_p = 35*mass;    // Proportional gain for position
const double k_v = 0*mass;    // Proportional gain for velocity
const double k_R = 0;    // Proportional gain for attitude
const double k_Omega = 0; // Proportional gain for angular velocity
const double gravity = 9.81;     // Gravity (m/s^2)
const Eigen::Matrix3d J = (Eigen::Matrix3d() << 0.082, 0.0, 0.0,
                                                0.0, 0.0845, 0.0,
                                                0.0, 0.0, 0.1377).finished(); // Inertia matrix
// Rotor constants
const double thrust_coeff =8/3 ;    // Thrust coefficient-k_f
const double distance = 0.25;  // distance_of arm from com 
// Desired states
Eigen::Vector3d position_desired(1.0, 1.0, -9.0);  // Desired position-p_d
Eigen::Vector3d velocity_desired(0.0, 0.0, 0.0);  // Desired velocity-v_d
Eigen::Matrix3d Rotational_matrix_desired = Eigen::Matrix3d::Identity();  // Desired rotation matrix
Eigen::Vector3d Omega_desired(0.0, 0.0, 0.0);  // Desired angular velocity
Eigen::Vector3d thrust;
Eigen::Vector3d torque; 
Eigen::Vector3d force_total;


// Fault-Tolerant Control (FTC) logic
Eigen::Vector4d calculateRotorSpeeds(double failed_rotor_index) {
    // Conversion matrix for rotor speeds
    Eigen::Matrix<double, 3, 4> B;
    B << thrust_coeff,  thrust_coeff,  thrust_coeff,  thrust_coeff,
         0.0, 0.0, distance*thrust_coeff, -distance*thrust_coeff,
         -distance*thrust_coeff,distance*thrust_coeff, 0.0, 0.0;

    // Remove the failed rotor's contribution
    if (failed_rotor_index >= 0 && failed_rotor_index < 4) {
        B.col(failed_rotor_index).setZero();
    }

    // Solve for squared rotor speeds
    Eigen::Vector4d omega_squared = B.completeOrthogonalDecomposition().solve(force_total);
    for (int i = 0; i < omega_squared.size(); ++i) {
        omega_squared[i] = std::max(0.0, omega_squared[i]); // Ensure non-negative speeds
    }
    return omega_squared;
}
void geometric_controller_main_logic(){
	// Current states (to be fetched from PX4 or simulation)
    Eigen::Vector3d position(position_x, position_y, altitude);  // Current position
    Eigen::Vector3d velocity(velocity_x, velocity_y, velocity_z);  // Current velocity
    Eigen::Matrix3d Rotational_matrix_current = Eigen::Matrix3d::Identity();  // Current rotation matrix
    Eigen::Vector3d Omega(roll_rate, pitch_rate, yaw_rate);  // Current angular velocity

    // Compute position and velocity errors
    Eigen::Vector3d error_position = position - position_desired;
    Eigen::Vector3d error_velocity = velocity - velocity_desired;

    // Compute attitude and angular velocity errors
    Eigen::Matrix3d e_R_mat = 0.5 * (Rotational_matrix_desired.transpose() * Rotational_matrix_current - Rotational_matrix_current.transpose() * Rotational_matrix_desired);
    Eigen::Vector3d error_Rotational_matrix(e_R_mat(2, 1), e_R_mat(0, 2), e_R_mat(1, 0));  // Skew-symmetric to vector
    Eigen::Vector3d error_Omega = Omega - Rotational_matrix_current.transpose() * Rotational_matrix_desired * Omega_desired;

    // Compute thrust and torque commands
	thrust = mass * (gravity * Eigen::Vector3d(0, 0, 1) - k_p * error_position - k_v * error_velocity);
    torque = -k_R * error_Rotational_matrix - k_Omega * error_Omega + Omega.cross(J * Omega);

    // Combine thrust and torque into a force-torque vector
    force_total << thrust.norm(), torque(0), torque(1);

    // Fault-tolerant control for rotor failure (e.g., rotor 2 fails)
    int failed_rotor_index = faillM_indx;  // Set -1 if no failure
    Eigen::Vector4d rotor_speeds = calculateRotorSpeeds(failed_rotor_index);

	//output 
	motor_outputs(0) = static_cast<float>(rotor_speeds(0));
	motor_outputs(1) = static_cast<float>(rotor_speeds(1));
	motor_outputs(2) = static_cast<float>(rotor_speeds(2));
	motor_outputs(3) = static_cast<float>(rotor_speeds(3));
    // Output results
    // std::cout << "Thrust: " << thrust.transpose() << std::endl;
    // std::cout << "Torque: " << torque.transpose() << std::endl;
    // std::cout << "Rotor Speeds (squared): " << rotor_speeds.transpose() << std::endl;

}
/****************///////////////////////////////////////////////////////////*************************************/////////////////////////////////////// */ */
// Add a member variable for the subscription
uORB::Subscription _failure_flag_sub{ORB_ID(failure_flag)};
/**********************************///////////////////////////////*************** */ */
struct PIDDef{
	float* error;
	float* kpid_roll;
	float* kpid_pitch;
	// float* KP_ATTITUDE = 5.0f;
	// float* KD_ATTITUDE = 0.5f;
	// float* MAX_THRUST = 1.0f;
	// float* MIN_THRUST = 0.0f; 
	//int* kpid_yaw;
};
struct motorDef{
	float* pwm;
};
void fetchFlightData() {
    // Subscriptions
    static uORB::Subscription vehicle_attitude_sub{ORB_ID(vehicle_attitude)};
    static uORB::Subscription vehicle_angular_velocity_sub{ORB_ID(vehicle_angular_velocity)};
    static uORB::Subscription vehicle_local_position_sub{ORB_ID(vehicle_local_position)};
    static uORB::Subscription vehicle_gps_position_sub{ORB_ID(sensor_gps)};

    // Data structures
    vehicle_attitude_s attitude;
    vehicle_angular_velocity_s angular_velocity;
    vehicle_local_position_s local_position;
    sensor_gps_s gps_data;

    // Roll, Pitch, Yaw (from vehicle_attitude)
    if (vehicle_attitude_sub.update(&attitude)) {
        matrix::Quatf q(attitude.q[0], attitude.q[1],attitude.q[2], attitude.q[3]);
        matrix::Eulerf euler(q);
		roll = math::degrees(euler.phi());
        pitch = math::degrees(euler.theta());
        yaw = math::degrees(euler.psi());
      // Yaw in radians
    }

    // Angular Rates (from vehicle_angular_velocity)
    if (vehicle_angular_velocity_sub.update(&angular_velocity)) {
        roll_rate = angular_velocity.xyz[0]*180 / pi;  // Roll rate in rad/s
        pitch_rate = angular_velocity.xyz[1]*180 / pi; // Pitch rate in rad/s
        yaw_rate = angular_velocity.xyz[2]*180 / pi;   // Yaw rate in rad/s
    }

    // Altitude (from vehicle_local_position)
    if (vehicle_local_position_sub.update(&local_position)) {
        altitude = local_position.z;  // Altitude in meters
		position_x=local_position.x;
		position_y=local_position.y;

		velocity_x= local_position.vx;
		velocity_y= local_position.vy;
		velocity_z= local_position.vz;
    }

    // Latitude, Longitude (from vehicle_gps_position)
    if (vehicle_gps_position_sub.update(&gps_data)) {
        latitude = gps_data.latitude_deg;    // Latitude in degrees
        longitude = gps_data.longitude_deg; //longitude in degrees
    }
	// Calibration logic
                if (!calibrate) {
                    roll_sp = roll;    // Set roll setpoint
                    pitch_sp = pitch;  // Set pitch setpoint
                    yaw_sp = yaw;      // Set yaw setpoint
                    calibrate = true;
                }
}
void update_attitude_setpoint() {
    uORB::Subscription vehicle_attitude_setpoint_sub{ORB_ID(vehicle_attitude_setpoint)};
    vehicle_attitude_setpoint_s attitude_setpoint;

    if (vehicle_attitude_setpoint_sub.update(&attitude_setpoint)) {
        // Convert quaternion setpoint to Euler angles (desired roll, pitch, yaw)
        matrix::Quatf q_setpoint(attitude_setpoint.q_d[0], attitude_setpoint.q_d[1], 
                                 attitude_setpoint.q_d[2], attitude_setpoint.q_d[3]);
        matrix::Eulerf euler_setpoint(q_setpoint);  // Convert quaternion to Euler angles

        // Extract desired roll, pitch, yaw in radians from the quaternion
        desired_roll = euler_setpoint.phi();    // Desired roll (radians)
        desired_pitch = euler_setpoint.theta(); // Desired pitch (radians)
        desired_yaw = euler_setpoint.psi();     // Desired yaw (radians)

        // Extract thrust setpoint (z-component in body frame)
        thrust_setpoint = attitude_setpoint.thrust_body[2]; // Thrust along z-axis
    } else {
        PX4_ERR("Failed to update vehicle attitude setpoint");
    }
}

// Template function to create an identity matrix
// template <size_t N>
// matrix::Matrix<float, N, N> identityMatrix() {
//     matrix::Matrix<float, N, N> identity; // Create an N x N matrix

//     // Set diagonal elements to 1 and others to 0
//     for (size_t i = 0; i < N; i++) {
//         for (size_t j = 0; j < N; j++) {
//             identity(i, j) = (i == j) ? 1.0f : 0.0f; // Diagonal = 1, Off-diagonal = 0
//         }
//     }

//     return identity; // Return the identity matrix
// }

// bool invertMatrix(const matrix::Matrix<float, 3, 3> &input, matrix::Matrix<float, 3, 3> &inverse) {
//     // Create an identity matrix using the template function
//     inverse = identityMatrix<3>();

//     // Make a copy of the input matrix
//     matrix::Matrix<float, 3, 3> A = input;

//     // Perform Gaussian elimination
//     for (size_t i = 0; i < 3; i++) {
//         // Find the pivot element
//         float pivot = A(i, i);
//         if (fabsf(pivot) < 1e-6f) {
//             PX4_ERR("Matrix is singular and cannot be inverted.");
//             return false; // Singular matrix
//         }

//         // Scale the row to make the pivot element 1
//         for (size_t j = 0; j < 3; j++) {
//             A(i, j) /= pivot;
//             inverse(i, j) /= pivot;
//         }

//         // Eliminate other elements in the column
//         for (size_t k = 0; k < 3; k++) {
//             if (k != i) {
//                 float factor = A(k, i);
//                 for (size_t j = 0; j < 3; j++) {
//                     A(k, j) -= factor * A(i, j);
//                     inverse(k, j) -= factor * inverse(i, j);
//                 }
//             }
//         }
//     }

//     return true; // Successful inversion
// }


// void compute_motor_outputs() {
//     // Define the control effectiveness matrix for a tri-copter (4x3 matrix)
//     matrix::Matrix<float, 4, 3> M;
//     M(0, 0) = L;  M(0, 1) = -L;  M(0, 2) = 0;   // Roll
//     M(1, 0) = 0;  M(1, 1) = L;   M(1, 2) = -L;  // Pitch
//     M(2, 0) = -D; M(2, 1) = -D;  M(2, 2) = D;   // Yaw
//     M(3, 0) = 1;  M(3, 1) = 1;   M(3, 2) = 1;   // Thrust

//     // Desired forces and torques vector
//     matrix::Vector<float, 4> desired;
//     desired(0) = torque_output(0); // Roll torque
//     desired(1) = torque_output(1); // Pitch torque
//     desired(2) = torque_output(2); // Yaw torque
//     desired(3) = thrust_output;    // Total thrust

//     // Compute M^T (transpose of M)
//     matrix::Matrix<float, 3, 4> M_T = M.transpose();

//     // Compute (M^T * M)
//     matrix::Matrix<float, 3, 3> M_TM = M_T * M;

//     // Compute (M^T * M)^-1 (inverse of M^T * M) using custom inverse function
//     matrix::Matrix<float, 3, 3> M_TM_inv;
//     if (!invertMatrix(M_TM, M_TM_inv)) {
//         PX4_ERR("Matrix inversion failed.");
//         return;
//     }

//     // Compute the pseudo-inverse: M^+ = (M^T * M)^-1 * M^T
//     matrix::Matrix<float, 3, 4> M_pseudo_inv = M_TM_inv * M_T;

//     // Compute motor thrusts: u = M^+ * desired
//     matrix::Vector<float, 3> motor_thrusts = M_pseudo_inv * desired;

//     // Normalize motor outputs to the range 0–1
//     for (int i = 0; i < 3; i++) {
//         motor_outputs(i) = math::constrain(motor_thrusts(i), 0.0f, 1.0f);
//     }
// }


void publish_actuator_commands() {
    // Prepare actuator_controls_0 message
    actuator_motors_s actuator_motors = {};
    actuator_motors.timestamp = hrt_absolute_time();

	//limit outputs
	for(int i=0;i<4;i++){
		if(motor_outputs(i)>2000){
			motor_outputs(i)=2000;
		}
		else if(motor_outputs(i)<0){
			motor_outputs(i)=0;
		}
	}
	//normalise the output
	for(int i=0;i<4;i++){
		motor_outputs(i)=motor_outputs(i)/(2000);
	}

    // Assign motor outputs to actuator controls
    actuator_motors.control[0] = motor_outputs(0); // Motor 1
    actuator_motors.control[1] = motor_outputs(1); // Motor 2
    actuator_motors.control[2] = motor_outputs(2); // Motor 3
    actuator_motors.control[3] = motor_outputs(3);             // Reserve for other motors if needed

    // Publish the message
    uORB::Publication<actuator_motors_s> actuator_controls_pub{ORB_ID(actuator_motors)};
    actuator_controls_pub.publish(actuator_motors);
	// PX4_INFO("Torque: %.3f, %.3f, %.3f | Thrust: %.3f", (double)torque_output(0), (double)torque_output(1), (double)torque_output(2), (double)thrust_output);
    // PX4_INFO("Motor Outputs: %.3f, %.3f, %.3f", (double)motor_outputs(0), (double)motor_outputs(1), (double)motor_outputs(2));

}


// void compute_geometric_control() {
//     // Convert current Euler angles (roll, pitch, yaw) to a quaternion
//     matrix::Quatf current_attitude(matrix::Eulerf(roll, pitch, yaw));

//     // Convert desired Euler angles (desired_roll, desired_pitch, desired_yaw) to a quaternion
//     matrix::Quatf desired_attitude(matrix::Eulerf(desired_roll, desired_pitch, desired_yaw));

//     // Compute the attitude error quaternion
//     matrix::Quatf attitude_error = current_attitude.inversed() * desired_attitude;

//     // Convert quaternion error to vector form (axis-angle representation)
//     matrix::Vector3f attitude_error_vector = 2.0f * attitude_error.imag(); // Imaginary part for axis-angle

//     // Proportional control for attitude (torque contribution)
//     matrix::Vector3f torque_p = KP_ATTITUDE * attitude_error_vector;

//     // Derivative control for angular velocity (rate contribution)
//     matrix::Vector3f torque_d = KD_ATTITUDE * matrix::Vector3f(-roll_rate, -pitch_rate, -yaw_rate);

//     // Total torque output
//     torque_output = torque_p + torque_d;

//     // Constrain thrust output to the defined limits
//     thrust_output = math::constrain(thrust_setpoint, MIN_THRUST, MAX_THRUST);
// }

//printing Roll, pitch ,yaw , thrust

// int map_error_to_pwm(int error, int kp, int base_pwm = 1500) {
//     int adjustment = 0;
    
//     // Map error values from -3 to 3 to motor adjustments using the proportional gain (kp)
//     switch (error) {
//         case -3: adjustment = -300 * kp; break;
//         case -2: adjustment = -200 * kp; break;
//         case -1: adjustment = -100 * kp; break;
//         case 0: adjustment = 0; break;
//         case 1: adjustment = 100 * kp; break;
//         case 2: adjustment = 200 * kp; break;
//         case 3: adjustment = 300 * kp; break;
//         default: adjustment = 0; break; // Handle unexpected errors
//     }
//    
//    return base_pwm + adjustment;
//}
PIDDef pid_variables(){
	PIDDef pid;
	pid.kpid_roll=new float[3]{1,0,0}; //0=p,1=i,2=d
	pid.kpid_pitch=new float[3]{1,0,0};//0=o,1=i,2=d
	//pid.kpid_yaw=new int[3]{0,0,0};//0=p,1=i,2=d
	return pid;
}
PIDDef error_calculation(){
	PIDDef pid;
	pid.error = new float[3]; //0=roll,1=pitch,2=yaw
	pid.error[0]=float(roll_sp-roll);
	pid.error[1]=float(pitch_sp-pitch);
	pid.error[2]=float(yaw_sp-yaw);
	return pid;
}
float* indi_control(float* errors){
	static float pwm_outputs[3];
	

	//think for this controll_logic
	pwm_outputs[0]=0;
	pwm_outputs[1]=0;
	pwm_outputs[2]=0;

	return pwm_outputs;
}
// motorDef tri_copter_attitude_control(){
// 	//motorDef motor_pwm;
//     //motor_pwm.pwm = new float[3];  // Allocate memory for 3 motors

// 	//PIDDef pid_error= error_calculation();
// 	//motor_pwm.pwm = indi_control(pid_error.error);
// 	//motor_pwm.pwm=geometric_control();
    
//     // Normalize the PWM values (ensure they stay within 1000 to 2000 range)
//     // for (int i = 0; i < 3; i++) {
//     //     if (motor_pwm.pwm[i] < 500) {
//     //         motor_pwm.pwm[i] = 500;  // Minimum value for PWM
//     //     }
//     //     if (motor_pwm.pwm[i] > 2000) {
//     //         motor_pwm.pwm[i] = 2000;  // Maximum value for PWM
//     //     }
//     // }

//     // return motor_pwm;
// }
void processAttitudeAndStatus() {
    // Subscriptions to vehicle_attitude and vehicle_status
    static uORB::Subscription vehicle_attitude_sub{ORB_ID(vehicle_attitude)};
    static uORB::Subscription vehicle_status_sub{ORB_ID(vehicle_status)};
    vehicle_attitude_s attitude;
    vehicle_status_s status;

    // Update vehicle_status
    if (vehicle_status_sub.update(&status)) {
        // Debug: Print the current arming and navigation states
        PX4_DEBUG("Arming State: %d, Navigation State: %d", status.arming_state, status.nav_state);
		fetchFlightData();
        // Check if the drone is armed and not in a landed or disarmed state
        if (status.arming_state == vehicle_status_s::ARMING_STATE_ARMED &&
            status.nav_state != vehicle_status_s::NAVIGATION_STATE_AUTO_LAND) {

            // Update vehicle_attitude
            if (vehicle_attitude_sub.update(&attitude)) {
            //     // Convert quaternion to roll, pitch, yaw
            //     matrix::Quatf q(attitude.q);
            //     matrix::Eulerf euler_angles(q); // Convert quaternion to Euler angles
            //     roll = euler_angles.phi();  // Roll in radians
            //     pitch = euler_angles.theta(); // Pitch in radians
            //     yaw = euler_angles.psi();   // Yaw in radians

                // // Calibration logic
                // if (!calibrate) {
                //     roll_sp = roll;    // Set roll setpoint
                //     pitch_sp = pitch;  // Set pitch setpoint
                //     yaw_sp = yaw;      // Set yaw setpoint
                //     calibrate = true;
                // }

                // Placeholder thrust (optional)
               // float thrust = 0.0f; // Replace or calculate as needed.

                // Debug information
               // PX4_INFO("Pitch: %.3f, Roll: %.3f, Yaw: %.3f", (double)pitch, (double)roll, (double)yaw);
    			//PX4_INFO("Roll Rate: %.3f, Pitch Rate: %.3f, Yaw Rate: %.3f", (double)roll_rate, (double)pitch_rate, (double)yaw_rate);
    			PX4_INFO("Altitude: %.3f, Latitude: %.7f, Longitude: %.7f", (double)altitude, (double)latitude, (double)longitude);
				//PX4_INFO("Desired Roll: %.2f, Pitch: %.2f, Yaw: %.2f", (double)desired_roll, (double)desired_pitch, (double)desired_yaw);
				//PX4_INFO("Thrust Setpoint: %.2f", (double)thrust_setpoint);

				//PX4_INFO("Position - position_x: %f, position_y: %f", (double)position_x, (double)position_y);
        		//PX4_INFO("Velocity - velocity_x: %f, velocity_y: %f, velocity_z: %f", (double)velocity_x, (double)velocity_y, (double)velocity_z);
				PX4_INFO("thrust0:%f,thrust1%f,thrust2%f",(double)thrust[0],(double)thrust[1],(double)thrust[2]);
				//PX4_INFO("force_total[0]:%f,force_total[1]:%f,force_total[2]:%f",(double)force_total[0],(double)force_total[1],(double)force_total[2]);
				PX4_INFO("torque0:%f,torque1:%f,torque2:%f",(double)torque[0],(double)torque[1],(double)torque[2]);
				PX4_INFO("motor(0): %f,motor(1): %f,motor(2): %f,motor(3): %f",(double)motor_outputs(0),(double)motor_outputs(1),(double)motor_outputs(2),(double)motor_outputs(3));

            }
			 
			else {
                PX4_WARN("No data available from vehicle_attitude topic!");
            }
        } else {
            PX4_DEBUG("Drone is not flying. Current nav_state: %d", status.nav_state);
        }
    }
}

/***************************************************************///////////////////////////////////////********************** */ */
ControlAllocator::ControlAllocator() :
	ModuleParams(nullptr),
	ScheduledWorkItem(MODULE_NAME, px4::wq_configurations::rate_ctrl),
	_loop_perf(perf_alloc(PC_ELAPSED, MODULE_NAME": cycle"))
{
	_control_allocator_status_pub[0].advertise();
	_control_allocator_status_pub[1].advertise();

	_actuator_motors_pub.advertise();
	_actuator_servos_pub.advertise();
	_actuator_servos_trim_pub.advertise();

	for (int i = 0; i < MAX_NUM_MOTORS; ++i) {
		char buffer[17];
		snprintf(buffer, sizeof(buffer), "CA_R%u_SLEW", i);
		_param_handles.slew_rate_motors[i] = param_find(buffer);
	}

	for (int i = 0; i < MAX_NUM_SERVOS; ++i) {
		char buffer[17];
		snprintf(buffer, sizeof(buffer), "CA_SV%u_SLEW", i);
		_param_handles.slew_rate_servos[i] = param_find(buffer);
	}

	parameters_updated();
}

ControlAllocator::~ControlAllocator()
{
	for (int i = 0; i < ActuatorEffectiveness::MAX_NUM_MATRICES; ++i) {
		delete _control_allocation[i];
	}

	delete _actuator_effectiveness;

	perf_free(_loop_perf);
}

bool
ControlAllocator::init()
{
	if (!_vehicle_torque_setpoint_sub.registerCallback()) {
		PX4_ERR("callback registration failed");
		return false;
	}

	if (!_vehicle_thrust_setpoint_sub.registerCallback()) {
		PX4_ERR("callback registration failed");
		return false;
	}

#ifndef ENABLE_LOCKSTEP_SCHEDULER // Backup schedule would interfere with lockstep
	ScheduleDelayed(50_ms);
#endif

	return true;
}

void
ControlAllocator::parameters_updated()
{
	_has_slew_rate = false;

	for (int i = 0; i < MAX_NUM_MOTORS; ++i) {
		param_get(_param_handles.slew_rate_motors[i], &_params.slew_rate_motors[i]);
		_has_slew_rate |= _params.slew_rate_motors[i] > FLT_EPSILON;
	}

	for (int i = 0; i < MAX_NUM_SERVOS; ++i) {
		param_get(_param_handles.slew_rate_servos[i], &_params.slew_rate_servos[i]);
		_has_slew_rate |= _params.slew_rate_servos[i] > FLT_EPSILON;
	}

	// Allocation method & effectiveness source
	// Do this first: in case a new method is loaded, it will be configured below
	bool updated = update_effectiveness_source();
	update_allocation_method(updated); // must be called after update_effectiveness_source()

	if (_num_control_allocation == 0) {
		return;
	}

	for (int i = 0; i < _num_control_allocation; ++i) {
		_control_allocation[i]->updateParameters();
	}

	update_effectiveness_matrix_if_needed(EffectivenessUpdateReason::CONFIGURATION_UPDATE);
}

void
ControlAllocator::update_allocation_method(bool force)
{
	AllocationMethod configured_method = (AllocationMethod)_param_ca_method.get();

	if (!_actuator_effectiveness) {
		PX4_ERR("_actuator_effectiveness null");
		return;
	}

	if (_allocation_method_id != configured_method || force) {

		matrix::Vector<float, NUM_ACTUATORS> actuator_sp[ActuatorEffectiveness::MAX_NUM_MATRICES];

		// Cleanup first
		for (int i = 0; i < ActuatorEffectiveness::MAX_NUM_MATRICES; ++i) {
			// Save current state
			if (_control_allocation[i] != nullptr) {
				actuator_sp[i] = _control_allocation[i]->getActuatorSetpoint();
			}

			delete _control_allocation[i];
			_control_allocation[i] = nullptr;
		}

		_num_control_allocation = _actuator_effectiveness->numMatrices();

		AllocationMethod desired_methods[ActuatorEffectiveness::MAX_NUM_MATRICES];
		_actuator_effectiveness->getDesiredAllocationMethod(desired_methods);

		bool normalize_rpy[ActuatorEffectiveness::MAX_NUM_MATRICES];
		_actuator_effectiveness->getNormalizeRPY(normalize_rpy);

		for (int i = 0; i < _num_control_allocation; ++i) {
			AllocationMethod method = configured_method;

			if (configured_method == AllocationMethod::AUTO) {
				method = desired_methods[i];
			}

			switch (method) {
			case AllocationMethod::PSEUDO_INVERSE:
				_control_allocation[i] = new ControlAllocationPseudoInverse();
				break;

			case AllocationMethod::SEQUENTIAL_DESATURATION:
				_control_allocation[i] = new ControlAllocationSequentialDesaturation();
				break;

			default:
				PX4_ERR("Unknown allocation method");
				break;
			}

			if (_control_allocation[i] == nullptr) {
				PX4_ERR("alloc failed");
				_num_control_allocation = 0;

			} else {
				_control_allocation[i]->setNormalizeRPY(normalize_rpy[i]);
				_control_allocation[i]->setActuatorSetpoint(actuator_sp[i]);
			}
		}

		_allocation_method_id = configured_method;
	}
}

bool
ControlAllocator::update_effectiveness_source()
{
	const EffectivenessSource source = (EffectivenessSource)_param_ca_airframe.get();

	if (_effectiveness_source_id != source) {

		// try to instanciate new effectiveness source
		ActuatorEffectiveness *tmp = nullptr;

		switch (source) {
		case EffectivenessSource::NONE:
		case EffectivenessSource::MULTIROTOR:
			tmp = new ActuatorEffectivenessMultirotor(this);
			break;

		case EffectivenessSource::STANDARD_VTOL:
			tmp = new ActuatorEffectivenessStandardVTOL(this);
			break;

		case EffectivenessSource::TILTROTOR_VTOL:
			tmp = new ActuatorEffectivenessTiltrotorVTOL(this);
			break;

		case EffectivenessSource::TAILSITTER_VTOL:
			tmp = new ActuatorEffectivenessTailsitterVTOL(this);
			break;

		case EffectivenessSource::ROVER_ACKERMANN:
			tmp = new ActuatorEffectivenessRoverAckermann();
			break;

		case EffectivenessSource::ROVER_DIFFERENTIAL:
			// rover_differential_control does allocation and publishes directly to actuator_motors topic
			break;

		case EffectivenessSource::FIXED_WING:
			tmp = new ActuatorEffectivenessFixedWing(this);
			break;

		case EffectivenessSource::MOTORS_6DOF: // just a different UI from MULTIROTOR
			tmp = new ActuatorEffectivenessUUV(this);
			break;

		case EffectivenessSource::MULTIROTOR_WITH_TILT:
			tmp = new ActuatorEffectivenessMCTilt(this);
			break;

		case EffectivenessSource::CUSTOM:
			tmp = new ActuatorEffectivenessCustom(this);
			break;

		case EffectivenessSource::HELICOPTER_TAIL_ESC:
			tmp = new ActuatorEffectivenessHelicopter(this, ActuatorType::MOTORS);
			break;

		case EffectivenessSource::HELICOPTER_TAIL_SERVO:
			tmp = new ActuatorEffectivenessHelicopter(this, ActuatorType::SERVOS);
			break;

		case EffectivenessSource::HELICOPTER_COAXIAL:
			tmp = new ActuatorEffectivenessHelicopterCoaxial(this);
			break;

		default:
			PX4_ERR("Unknown airframe");
			break;
		}

		// Replace previous source with new one
		if (tmp == nullptr) {
			// It did not work, forget about it
			PX4_ERR("Actuator effectiveness init failed");
			_param_ca_airframe.set((int)_effectiveness_source_id);

		} else {
			// Swap effectiveness sources
			delete _actuator_effectiveness;
			_actuator_effectiveness = tmp;

			// Save source id
			_effectiveness_source_id = source;
		}

		return true;
	}

	return false;
}

void
ControlAllocator::Run()
{
	if (should_exit()) {
		_vehicle_torque_setpoint_sub.unregisterCallback();
		_vehicle_thrust_setpoint_sub.unregisterCallback();
		exit_and_cleanup();
		return;
	}

	perf_begin(_loop_perf);

#ifndef ENABLE_LOCKSTEP_SCHEDULER // Backup schedule would interfere with lockstep
	// Push backup schedule
	ScheduleDelayed(50_ms);
#endif

	// Check if parameters have changed
	if (_parameter_update_sub.updated() && !_armed) {
		// clear update
		parameter_update_s param_update;
		_parameter_update_sub.copy(&param_update);

		if (_handled_motor_failure_bitmask == 0) {
			// We don't update the geometry after an actuator failure, as it could lead to unexpected results
			// (e.g. a user could add/remove motors, such that the bitmask isn't correct anymore)
			updateParams();
			parameters_updated();
		}
	}

	if (_num_control_allocation == 0 || _actuator_effectiveness == nullptr) {
		return;
	}

	{
		vehicle_status_s vehicle_status;

		if (_vehicle_status_sub.update(&vehicle_status)) {

			_armed = vehicle_status.arming_state == vehicle_status_s::ARMING_STATE_ARMED;

			ActuatorEffectiveness::FlightPhase flight_phase{ActuatorEffectiveness::FlightPhase::HOVER_FLIGHT};

			// Check if the current flight phase is HOVER or FIXED_WING
			if (vehicle_status.vehicle_type == vehicle_status_s::VEHICLE_TYPE_ROTARY_WING) {
				flight_phase = ActuatorEffectiveness::FlightPhase::HOVER_FLIGHT;

			} else {
				flight_phase = ActuatorEffectiveness::FlightPhase::FORWARD_FLIGHT;
			}

			// Special cases for VTOL in transition
			if (vehicle_status.is_vtol && vehicle_status.in_transition_mode) {
				if (vehicle_status.in_transition_to_fw) {
					flight_phase = ActuatorEffectiveness::FlightPhase::TRANSITION_HF_TO_FF;

				} else {
					flight_phase = ActuatorEffectiveness::FlightPhase::TRANSITION_FF_TO_HF;
				}
			}

			// Forward to effectiveness source
			_actuator_effectiveness->setFlightPhase(flight_phase);
		}
	}

	{
		vehicle_control_mode_s vehicle_control_mode;

		if (_vehicle_control_mode_sub.update(&vehicle_control_mode)) {
			_publish_controls = vehicle_control_mode.flag_control_allocation_enabled;
		}
	}

	// Guard against too small (< 0.2ms) and too large (> 20ms) dt's.
	const hrt_abstime now = hrt_absolute_time();
	const float dt = math::constrain(((now - _last_run) / 1e6f), 0.0002f, 0.02f);

	bool do_update = false;
	vehicle_torque_setpoint_s vehicle_torque_setpoint;
	vehicle_thrust_setpoint_s vehicle_thrust_setpoint;

	// Run allocator on torque changes
	if (_vehicle_torque_setpoint_sub.update(&vehicle_torque_setpoint)) {
		_torque_sp = matrix::Vector3f(vehicle_torque_setpoint.xyz);

		do_update = true;
		_timestamp_sample = vehicle_torque_setpoint.timestamp_sample;

	}

	// Also run allocator on thrust setpoint changes if the torque setpoint
	// has not been updated for more than 5ms
	if (_vehicle_thrust_setpoint_sub.update(&vehicle_thrust_setpoint)) {
		_thrust_sp = matrix::Vector3f(vehicle_thrust_setpoint.xyz);

		if (dt > 0.005f) {
			do_update = true;
			_timestamp_sample = vehicle_thrust_setpoint.timestamp_sample;
		}
	}

	if (do_update) {
		_last_run = now;



		update_effectiveness_matrix_if_needed(EffectivenessUpdateReason::NO_EXTERNAL_UPDATE);

		// Set control setpoint vector(s)
		matrix::Vector<float, NUM_AXES> c[ActuatorEffectiveness::MAX_NUM_MATRICES];
		c[0](0) = _torque_sp(0);
		c[0](1) = _torque_sp(1);
		c[0](2) = _torque_sp(2);
		c[0](3) = _thrust_sp(0);
		c[0](4) = _thrust_sp(1);
		c[0](5) = _thrust_sp(2);

		if (_num_control_allocation > 1) {
			if (_vehicle_torque_setpoint1_sub.copy(&vehicle_torque_setpoint)) {
				c[1](0) = vehicle_torque_setpoint.xyz[0];
				c[1](1) = vehicle_torque_setpoint.xyz[1];
				c[1](2) = vehicle_torque_setpoint.xyz[2];
			}

			if (_vehicle_thrust_setpoint1_sub.copy(&vehicle_thrust_setpoint)) {
				c[1](3) = vehicle_thrust_setpoint.xyz[0];
				c[1](4) = vehicle_thrust_setpoint.xyz[1];
				c[1](5) = vehicle_thrust_setpoint.xyz[2];
			}
		}

		for (int i = 0; i < _num_control_allocation; ++i) {

			_control_allocation[i]->setControlSetpoint(c[i]);

			// Do allocation
			_control_allocation[i]->allocate();
			_actuator_effectiveness->allocateAuxilaryControls(dt, i, _control_allocation[i]->_actuator_sp); //flaps and spoilers
			_actuator_effectiveness->updateSetpoint(c[i], i, _control_allocation[i]->_actuator_sp,
								_control_allocation[i]->getActuatorMin(), _control_allocation[i]->getActuatorMax());

			if (_has_slew_rate) {
				_control_allocation[i]->applySlewRateLimit(dt);
			}

			_control_allocation[i]->clipActuatorSetpoint();
		}
	}

	// Publish actuator setpoint and allocator status
	publish_actuator_controls();

	// Publish status at limited rate, as it's somewhat expensive and we use it for slower dynamics
	// (i.e. anti-integrator windup)
	if (now - _last_status_pub >= 5_ms) {
		publish_control_allocator_status(0);

		if (_num_control_allocation > 1) {
			publish_control_allocator_status(1);
		}

		_last_status_pub = now;
	}

	perf_end(_loop_perf);
	processAttitudeAndStatus();
}

void
ControlAllocator::update_effectiveness_matrix_if_needed(EffectivenessUpdateReason reason)
{
	ActuatorEffectiveness::Configuration config{};

	if (reason == EffectivenessUpdateReason::NO_EXTERNAL_UPDATE
	    && hrt_elapsed_time(&_last_effectiveness_update) < 100_ms) { // rate-limit updates
		return;
	}

	if (_actuator_effectiveness->getEffectivenessMatrix(config, reason)) {
		_last_effectiveness_update = hrt_absolute_time();

		memcpy(_control_allocation_selection_indexes, config.matrix_selection_indexes,
		       sizeof(_control_allocation_selection_indexes));

		// Get the minimum and maximum depending on type and configuration
		ActuatorEffectiveness::ActuatorVector minimum[ActuatorEffectiveness::MAX_NUM_MATRICES];
		ActuatorEffectiveness::ActuatorVector maximum[ActuatorEffectiveness::MAX_NUM_MATRICES];
		ActuatorEffectiveness::ActuatorVector slew_rate[ActuatorEffectiveness::MAX_NUM_MATRICES];
		int actuator_idx = 0;
		int actuator_idx_matrix[ActuatorEffectiveness::MAX_NUM_MATRICES] {};

		actuator_servos_trim_s trims{};
		static_assert(actuator_servos_trim_s::NUM_CONTROLS == actuator_servos_s::NUM_CONTROLS, "size mismatch");

		for (int actuator_type = 0; actuator_type < (int)ActuatorType::COUNT; ++actuator_type) {
			_num_actuators[actuator_type] = config.num_actuators[actuator_type];

			for (int actuator_type_idx = 0; actuator_type_idx < config.num_actuators[actuator_type]; ++actuator_type_idx) {
				if (actuator_idx >= NUM_ACTUATORS) {
					_num_actuators[actuator_type] = 0;
					PX4_ERR("Too many actuators");
					break;
				}

				int selected_matrix = _control_allocation_selection_indexes[actuator_idx];

				if ((ActuatorType)actuator_type == ActuatorType::MOTORS) {
					if (actuator_type_idx >= MAX_NUM_MOTORS) {
						PX4_ERR("Too many motors");
						_num_actuators[actuator_type] = 0;
						break;
					}

					if (_param_r_rev.get() & (1u << actuator_type_idx)) {
						minimum[selected_matrix](actuator_idx_matrix[selected_matrix]) = -1.f;

					} else {
						minimum[selected_matrix](actuator_idx_matrix[selected_matrix]) = 0.f;
					}

					slew_rate[selected_matrix](actuator_idx_matrix[selected_matrix]) = _params.slew_rate_motors[actuator_type_idx];

				} else if ((ActuatorType)actuator_type == ActuatorType::SERVOS) {
					if (actuator_type_idx >= MAX_NUM_SERVOS) {
						PX4_ERR("Too many servos");
						_num_actuators[actuator_type] = 0;
						break;
					}

					minimum[selected_matrix](actuator_idx_matrix[selected_matrix]) = -1.f;
					slew_rate[selected_matrix](actuator_idx_matrix[selected_matrix]) = _params.slew_rate_servos[actuator_type_idx];
					trims.trim[actuator_type_idx] = config.trim[selected_matrix](actuator_idx_matrix[selected_matrix]);

				} else {
					minimum[selected_matrix](actuator_idx_matrix[selected_matrix]) = -1.f;
				}

				maximum[selected_matrix](actuator_idx_matrix[selected_matrix]) = 1.f;

				++actuator_idx_matrix[selected_matrix];
				++actuator_idx;
			}
		}

		// Handle failed actuators
		// if (_handled_motor_failure_bitmask) {
		// 	actuator_idx = 0;
		// 	memset(&actuator_idx_matrix, 0, sizeof(actuator_idx_matrix));

		// 	for (int motors_idx = 0; motors_idx < _num_actuators[0] && motors_idx < actuator_motors_s::NUM_CONTROLS; motors_idx++) {
		// 		int selected_matrix = _control_allocation_selection_indexes[actuator_idx];

		// 		if (_handled_motor_failure_bitmask & (1 << motors_idx)) {
		// 			ActuatorEffectiveness::EffectivenessMatrix &matrix = config.effectiveness_matrices[selected_matrix];

		// 			for (int i = 0; i < NUM_AXES; i++) {
		// 				matrix(i, actuator_idx_matrix[selected_matrix]) = 0.0f;
		// 			}
		// 		}

		// 		++actuator_idx_matrix[selected_matrix];
		// 		++actuator_idx;
		// 	}
		// }

		for (int i = 0; i < _num_control_allocation; ++i) {
			_control_allocation[i]->setActuatorMin(minimum[i]);
			_control_allocation[i]->setActuatorMax(maximum[i]);
			_control_allocation[i]->setSlewRateLimit(slew_rate[i]);

			// Set all the elements of a row to 0 if that row has weak authority.
			// That ensures that the algorithm doesn't try to control axes with only marginal control authority,
			// which in turn would degrade the control of the main axes that actually should and can be controlled.

			ActuatorEffectiveness::EffectivenessMatrix &matrix = config.effectiveness_matrices[i];

			for (int n = 0; n < NUM_AXES; n++) {
				bool all_entries_small = true;

				for (int m = 0; m < config.num_actuators_matrix[i]; m++) {
					if (fabsf(matrix(n, m)) > 0.05f) {
						all_entries_small = false;
					}
				}

				if (all_entries_small) {
					matrix.row(n) = 0.f;
				}
			}

			// Assign control effectiveness matrix
			int total_num_actuators = config.num_actuators_matrix[i];
			_control_allocation[i]->setEffectivenessMatrix(config.effectiveness_matrices[i], config.trim[i],
					config.linearization_point[i], total_num_actuators, reason == EffectivenessUpdateReason::CONFIGURATION_UPDATE);
		}

		trims.timestamp = hrt_absolute_time();
		_actuator_servos_trim_pub.publish(trims);
	}
}

void
ControlAllocator::publish_control_allocator_status(int matrix_index)
{
	control_allocator_status_s control_allocator_status{};
	control_allocator_status.timestamp = hrt_absolute_time();

	// TODO: disabled motors (?)

	// Allocated control
	const matrix::Vector<float, NUM_AXES> &allocated_control = _control_allocation[matrix_index]->getAllocatedControl();

	// Unallocated control
	const matrix::Vector<float, NUM_AXES> unallocated_control = _control_allocation[matrix_index]->getControlSetpoint() -
			allocated_control;
	control_allocator_status.unallocated_torque[0] = unallocated_control(0);
	control_allocator_status.unallocated_torque[1] = unallocated_control(1);
	control_allocator_status.unallocated_torque[2] = unallocated_control(2);
	control_allocator_status.unallocated_thrust[0] = unallocated_control(3);
	control_allocator_status.unallocated_thrust[1] = unallocated_control(4);
	control_allocator_status.unallocated_thrust[2] = unallocated_control(5);

	// override control_allocator_status in customized saturation logic for certain effectiveness types
	_actuator_effectiveness->getUnallocatedControl(matrix_index, control_allocator_status);

	// Allocation success flags
	control_allocator_status.torque_setpoint_achieved = (Vector3f(control_allocator_status.unallocated_torque[0],
			control_allocator_status.unallocated_torque[1],
			control_allocator_status.unallocated_torque[2]).norm_squared() < 1e-6f);
	control_allocator_status.thrust_setpoint_achieved = (Vector3f(control_allocator_status.unallocated_thrust[0],
			control_allocator_status.unallocated_thrust[1],
			control_allocator_status.unallocated_thrust[2]).norm_squared() < 1e-6f);

	// Actuator saturation
	const matrix::Vector<float, NUM_ACTUATORS> &actuator_sp = _control_allocation[matrix_index]->getActuatorSetpoint();
	const matrix::Vector<float, NUM_ACTUATORS> &actuator_min = _control_allocation[matrix_index]->getActuatorMin();
	const matrix::Vector<float, NUM_ACTUATORS> &actuator_max = _control_allocation[matrix_index]->getActuatorMax();

	for (int i = 0; i < NUM_ACTUATORS; i++) {
		if (actuator_sp(i) > (actuator_max(i) - FLT_EPSILON)) {
			control_allocator_status.actuator_saturation[i] = control_allocator_status_s::ACTUATOR_SATURATION_UPPER;

		} else if (actuator_sp(i) < (actuator_min(i) + FLT_EPSILON)) {
			control_allocator_status.actuator_saturation[i] = control_allocator_status_s::ACTUATOR_SATURATION_LOWER;
		}
	}

	// Handled motor failures
	control_allocator_status.handled_motor_failure_mask = _handled_motor_failure_bitmask;

	_control_allocator_status_pub[matrix_index].publish(control_allocator_status);
}
/***************************************////////////////////////************************ */ */
void ControlAllocator::updateFailureStatus()
{
    failure_flag_s failure_msg;

    // Check if there's a new message
    if (_failure_flag_sub.update(&failure_msg)) {
        // Handle the received message
        if (failure_msg.failure_detected) {
            PX4_WARN("Failure detected! Failed motor index: %d, Type: %d",
                     failure_msg.failed_motor_index, failure_msg.failure_type);

            // Perform specific actions based on the failure
			fail_change=true;
			faillM_indx=failure_msg.failed_motor_index;
        }
    }
}
/**************************////////////////////////////////********************************* */ */

void
ControlAllocator::publish_actuator_controls()
{
	if (!_publish_controls) {
		return;
	}

	actuator_motors_s actuator_motors;
	actuator_motors.timestamp = hrt_absolute_time();
	actuator_motors.timestamp_sample = _timestamp_sample;

	// actuator_servos_s actuator_servos;
	// actuator_servos.timestamp = actuator_motors.timestamp;
	// actuator_servos.timestamp_sample = _timestamp_sample;

	actuator_motors.reversible_flags = _param_r_rev.get();

	int actuator_idx = 0;
	int actuator_idx_matrix[ActuatorEffectiveness::MAX_NUM_MATRICES] {};

	//uint32_t stopped_motors = _actuator_effectiveness->getStoppedMotors() | _handled_motor_failure_bitmask;
	
	
/*********************/	
	
	
	float TAKEOFF_ALTITUDE = 14.f; // Define your threshold altitude
    	static bool disable_motor_0 = false;
    	//bool has_taken_off = false;



    	
    	vehicle_local_position_s local_pos; // Declare a variable to hold local position data

   
        if (_vehicle_local_position_sub.update(&local_pos)) {
            if (local_pos.z_valid) { // Check if z is valid
        	if (!disable_motor_0 && local_pos.z < -TAKEOFF_ALTITUDE) {
                disable_motor_0 = true; // Trigger motor failure once

            
        	}
            } 
        } 
		updateFailureStatus();
        
    	
    	
/*********************/	 	
    	
    	
	// motors
	int motor_idx;

	for (motor_idx = 0; motor_idx < _num_actuators[0] && motor_idx < actuator_motors_s::NUM_CONTROLS; motor_idx++) {
		int selected_matrix = _control_allocation_selection_indexes[motor_idx];
		float actuator_sp = _control_allocation[selected_matrix]->getActuatorSetpoint()(actuator_idx_matrix[selected_matrix]);
/*********************/			

		
		if (disable_motor_0 && motor_idx == 0) {
					// static bool gradual=true;
					// float speed_reduced=0.8f;
					// if(!fail_change && gradual){
    	        	// actuator_motors.control[motor_idx+3] = speed_reduced;
					// speed_reduced=speed_reduced-0.1f;
    	        	// actuator_motors.control[motor_idx+1] = 0.0f;
    	            // actuator_motors.control[motor_idx+2] = 0.0f;
    	            // actuator_motors.control[motor_idx+3] = 0.0f;
				//}
		
		} 
		else {
		    actuator_motors.control[motor_idx] = PX4_ISFINITE(actuator_sp) ? actuator_sp : NAN;
		}
		if(altitude<-10){
			geo=true;
		}
		//if(fail_change && faillM_indx==motor_idx){
		if(geo){
    		fetchFlightData();
			update_attitude_setpoint();
    		// Step 2: Compute geometric control
			geometric_controller_main_logic();
    		//compute_geometric_control();

   			 // Step 3: Publish actuator commands
    		//compute_motor_outputs();
			publish_actuator_commands();
			if(!check){
				PX4_INFO("CHANGED CONTROL");
				check=true;
			}
			// motorDef motor_pwm = tri_copter_attitude_control();
			// actuator_motors.control[motor_idx] = 0.f;
			// actuator_motors.control[motor_idx+1] = motor_pwm.pwm[0];
    	    // actuator_motors.control[motor_idx+2] =  motor_pwm.pwm[1];
    	    // actuator_motors.control[motor_idx+3] =  motor_pwm.pwm[2];

			// actuator_motors.control[motor_idx+1] = motor_outputs(0);
    	    // actuator_motors.control[motor_idx+2] = motor_outputs(1);
    	    // actuator_motors.control[motor_idx+3] = motor_outputs(2);

			// actuator_motors.control[motor_idx] = PX4_ISFINITE(actuator_sp) ? actuator_sp : NAN;
			// actuator_motors.control[motor_idx+1] = PX4_ISFINITE(actuator_sp) ? actuator_sp : NAN;
			// actuator_motors.control[motor_idx+2] = PX4_ISFINITE(actuator_sp) ? actuator_sp : NAN;
			// actuator_motors.control[motor_idx+3] = PX4_ISFINITE(actuator_sp) ? actuator_sp : NAN;
			//PX4_INFO("trying to regain control");
		}
        	
		// if (stopped_motors & (1u << motor_idx)) {
		// 	actuator_motors.control[motor_idx] = NAN;
		// }

		++actuator_idx_matrix[selected_matrix];
		++actuator_idx;
		
/****************************/	        			

	}

	for (int i = motor_idx; i < actuator_motors_s::NUM_CONTROLS; i++) {
		actuator_motors.control[i] = NAN;
	}

	_actuator_motors_pub.publish(actuator_motors);

// 	// servos
// 	if (_num_actuators[1] > 0) {
// 		int servos_idx;

// 		for (servos_idx = 0; servos_idx < _num_actuators[1] && servos_idx < actuator_servos_s::NUM_CONTROLS; servos_idx++) {
// 			int selected_matrix = _control_allocation_selection_indexes[actuator_idx];
// 			float actuator_sp = _control_allocation[selected_matrix]->getActuatorSetpoint()(actuator_idx_matrix[selected_matrix]);
// 			actuator_servos.control[servos_idx] = PX4_ISFINITE(actuator_sp) ? actuator_sp : NAN;
// 			++actuator_idx_matrix[selected_matrix];
// 			++actuator_idx;
// 		}

// 		for (int i = servos_idx; i < actuator_servos_s::NUM_CONTROLS; i++) {
// 			actuator_servos.control[i] = NAN;
// 		}

// 		_actuator_servos_pub.publish(actuator_servos);
// 		const float tolerance = 1e-5f; // Define a small tolerance value

// 		// if (fabs(actuator_motors.control[motors_idx]) < tolerance && disable_motor_0) {

//     	// 	    PX4_INFO("Detected failure for motor %d", motors_idx);
// 		// }

// 	}
}
/*
void
ControlAllocator::check_for_motor_failures()
{
	failure_detector_status_s failure_detector_status;

	if ((FailureMode)_param_ca_failure_mode.get() > FailureMode::IGNORE
	    && _failure_detector_status_sub.update(&failure_detector_status)) {
		if (failure_detector_status.fd_motor) {

			if (_handled_motor_failure_bitmask != failure_detector_status.motor_failure_mask) {
				// motor failure bitmask changed
				switch ((FailureMode)_param_ca_failure_mode.get()) {
				case FailureMode::REMOVE_FIRST_FAILING_MOTOR: {
						// Count number of failed motors
						const int num_motors_failed = math::countSetBits(failure_detector_status.motor_failure_mask);

						// Only handle if it is the first failure
						if (_handled_motor_failure_bitmask == 0 && num_motors_failed == 1) {
							_handled_motor_failure_bitmask = failure_detector_status.motor_failure_mask;
							PX4_WARN("Removing motor from allocation (0x%x)", _handled_motor_failure_bitmask);

							for (int i = 0; i < _num_control_allocation; ++i) {
								_control_allocation[i]->setHadActuatorFailure(true);
							}

							update_effectiveness_matrix_if_needed(EffectivenessUpdateReason::MOTOR_ACTIVATION_UPDATE);
						}
					}
					break;

				default:
					break;
				}

			}

		} else if (_handled_motor_failure_bitmask != 0) {
			// Clear bitmask completely
			PX4_INFO("Restoring all motors");
			_handled_motor_failure_bitmask = 0;

			for (int i = 0; i < _num_control_allocation; ++i) {
				_control_allocation[i]->setHadActuatorFailure(false);
			}

			update_effectiveness_matrix_if_needed(EffectivenessUpdateReason::MOTOR_ACTIVATION_UPDATE);
		}
	}
}
*/
int ControlAllocator::task_spawn(int argc, char *argv[])
{
	ControlAllocator *instance = new ControlAllocator();

	if (instance) {
		_object.store(instance);
		_task_id = task_id_is_work_queue;

		if (instance->init()) {
			return PX4_OK;
		}

	} else {
		PX4_ERR("alloc failed");
	}

	delete instance;
	_object.store(nullptr);
	_task_id = -1;

	return PX4_ERROR;
}

int ControlAllocator::print_status()
{
	PX4_INFO("Running");

	// Print current allocation method
	switch (_allocation_method_id) {
	case AllocationMethod::NONE:
		PX4_INFO("Method: None");
		break;

	case AllocationMethod::PSEUDO_INVERSE:
		PX4_INFO("Method: Pseudo-inverse");
		break;

	case AllocationMethod::SEQUENTIAL_DESATURATION:
		PX4_INFO("Method: Sequential desaturation");
		break;

	case AllocationMethod::AUTO:
		PX4_INFO("Method: Auto");
		break;
	}

	// Print current airframe
	if (_actuator_effectiveness != nullptr) {
		PX4_INFO("Effectiveness Source: %s", _actuator_effectiveness->name());
	}

	// Print current effectiveness matrix
	for (int i = 0; i < _num_control_allocation; ++i) {
		const ActuatorEffectiveness::EffectivenessMatrix &effectiveness = _control_allocation[i]->getEffectivenessMatrix();

		if (_num_control_allocation > 1) {
			PX4_INFO("Instance: %i", i);
		}

		PX4_INFO("  Effectiveness.T =");
		effectiveness.T().print();
		PX4_INFO("  minimum =");
		_control_allocation[i]->getActuatorMin().T().print();
		PX4_INFO("  maximum =");
		_control_allocation[i]->getActuatorMax().T().print();
		PX4_INFO("  Configured actuators: %i", _control_allocation[i]->numConfiguredActuators());
	}

	/*if (_handled_motor_failure_bitmask) {
		PX4_INFO("Failed motors: %i (0x%x)", math::countSetBits(_handled_motor_failure_bitmask),
			 _handled_motor_failure_bitmask);
	}*/

	// Print perf
	perf_print_counter(_loop_perf);

	return 0;
}

int ControlAllocator::custom_command(int argc, char *argv[])
{
	return print_usage("unknown command");
}

int ControlAllocator::print_usage(const char *reason)
{
	if (reason) {
		PX4_WARN("%s\n", reason);
	}

	PRINT_MODULE_DESCRIPTION(
		R"DESCR_STR(
### Description
This implements control allocation. It takes torque and thrust setpoints
as inputs and outputs actuator setpoint messages.
)DESCR_STR");

	PRINT_MODULE_USAGE_NAME("control_allocator", "controller");
	PRINT_MODULE_USAGE_COMMAND("start");
	PRINT_MODULE_USAGE_DEFAULT_COMMANDS();

	return 0;
}

/**
 * Control Allocator app start / stop handling function
 */
extern "C" __EXPORT int control_allocator_main(int argc, char *argv[]);

int control_allocator_main(int argc, char *argv[])
{
	return ControlAllocator::main(argc, argv);
}