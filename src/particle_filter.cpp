/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
std::default_random_engine gen;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  // Setting the number of particles
  num_particles = 100;

  //@param std[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m], standard deviation of yaw [rad]]
  double std_x = std[0];
  double std_y = std[1];
  double std_theta = std[2];
  
  // Creating Normal/Gaussian Distribution for x, y and theta
  std::normal_distribution<double> normal_dist_x(x, std_x);
  std::normal_distribution<double> normal_dist_y(y, std_y);
  std::normal_distribution<double> normal_dist_theta(theta, std_theta);
  
  // Initializing all particles to first position and adding random Gaussian noise to each particle
  for (int i = 0; i < num_particles; ++i) {
    //Using the struct Particle (declared in particle_filter.h)
    Particle particle;
    particle.id = i;
    //@param x Initial x position [m] (simulated estimate from GPS); normal_dist_ is the normal distribution (and 'gen' is the random engine)
    particle.x = normal_dist_x(gen);
    //@param y Initial y position [m]
    particle.y = normal_dist_y(gen);
    //@param theta Initial orientation [rad]
    particle.theta = normal_dist_theta(gen);
    //Initializing all weights to 1
    particle.weight = 1.0;
    particles.push_back(particle);
  }
  
  //flag if filter is initialized
  is_initialized = true;
  
  //std::cout << "Debug 1 - Initialization!" << std::endl;
  
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  // Creating Normal/Gaussian distribution for x, y and theta noises
  //@param std_pos[] Array of dimension 3 [standard deviation of x [m], standard deviation of y [m], standard deviation of yaw [rad]]
  std::normal_distribution<double> normal_dist_x(0, std_pos[0]);
  std::normal_distribution<double> normal_dist_y(0, std_pos[1]);
  std::normal_distribution<double> normal_dist_theta(0, std_pos[2]);

  for (int i = 0; i < num_particles; i++) {
    //Depending on if yaw_rate is zero or not, we're gonna add different measurements here... to aboid dividing by zero, I'm gonna compare yaw to a number very closo to zero instead
    if (fabs(yaw_rate) < 0.0000001) {
      // motion model when yaw rate = 0
      particles[i].x += velocity * delta_t * cos(particles[i].theta);
      particles[i].y += velocity * delta_t * sin(particles[i].theta);
      //Because, in this case yaw_rate equalts to 0, theta will stay the same
    } else {
      // motion model when yaw rate != 0
      particles[i].x += (velocity/yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
      particles[i].y += (velocity/yaw_rate) * (cos(particles[i].theta)-cos(particles[i].theta + yaw_rate * delta_t));
      particles[i].theta += yaw_rate * delta_t;
    }                                     
      //Adding gaussian noise to movement to account for uncertainty 
      particles[i].x += normal_dist_x(gen);
      particles[i].y += normal_dist_y(gen);
      particles[i].theta += normal_dist_theta(gen);                                        
  }
  //std::cout << "Debug 2 - Prediction!" << std::endl;
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  for (unsigned int i = 0; i < observations.size(); i++) {
    
    // min distance to max
    double min_dist = std::numeric_limits<double>::max();
    
    int map_id = -1;
    
    for (unsigned int j = 0; j < predicted.size(); j++) {
      
      // get distance bteween current/predicted landmarks
      double current_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      
      //find predicted landmark nearest the current observed landmark
      if (current_dist < min_dist) {
        min_dist = current_dist;
        map_id = predicted[j].id;
      }
    }
 
    // set the observation's id to the nearest predicted landmark's id
    observations[i].id = map_id;
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  // Following the recommendation from the Udacity's instructor (the mentor Neha V), in this phase I'll basically build an algorithm that ensures the following demands:
  // TRANSFORM each observation marker from the vehicle's coordinates to the map's coordinates
  // ENSURE map landmarks are inside sensor range -> Here, I'll loop through the map landmarks and find the ones which are within sensor range
  // Nearest Neighbor Data Association -> Here, I'll do data asociation between "Step 1" and "Step 2" vectors
  // Compute WEIGHT of particle
  
  //Loop for each particle
  for(int i = 0; i < num_particles; i++) {
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;
  
    
    //Create a vector to hold the map landmark locations predicted to be within sensor range of the particle
    vector<LandmarkObs> predictions;
    
    //Loop for each landmark
    for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      //Get id and x,y coordinates
      float lm_x = map_landmarks.landmark_list[j].x_f;
      float lm_y = map_landmarks.landmark_list[j].y_f;
      int lm_id = map_landmarks.landmark_list[j].id_i;
    
      // Here we're only going to consider landmarks with values inside the sensor range of the particle
      if(fabs(lm_x - x) <= sensor_range && fabs(lm_y - y) <= sensor_range) {
        predictions.push_back(LandmarkObs{ lm_id, lm_x, lm_y });
      }
    }
    
    vector<LandmarkObs> trans_os;
    for (unsigned int j = 0; j < observations.size(); j++) {
      double tx = x + cos(theta) * observations[j].x - sin(theta) * observations[j].y;
      double ty = y + sin(theta) * observations[j].x + cos(theta) * observations[j].y;
      trans_os.push_back(LandmarkObs{observations[j].id, tx, ty});
    }
    
    //Applying Data Association for predictions and transformed observations
    dataAssociation(predictions, trans_os);
    particles[i].weight = 1.0;
    for (unsigned int j = 0; j < trans_os.size(); j++) {
      double observed_x, observed_y, predicted_x, predicted_y;
      observed_x = trans_os[j].x;
      observed_y = trans_os[j].y;
      int data_association_prediction = trans_os[j].id;
      
      for (unsigned int k = 0; k < predictions.size(); k++) {
        if (predictions[k].id == data_association_prediction) {
          predicted_x = predictions[k].x;
          predicted_y = predictions[k].y;
        }
      }
      
      //Using multivariate Gaussian on predictions
      double s_x = std_landmark[0];
      double s_y = std_landmark[1];
      double obs_w = ( 1/(2*M_PI*s_x*s_y)) * exp( -( pow(predicted_x-observed_x,2)/(2*pow(s_x, 2)) + (pow(predicted_y-observed_y,2)/(2*pow(s_y, 2))) ) );
      
      //Product of this obsevation weight with total observations weight
      particles[i].weight *= obs_w;
    }
  
  }
  
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  
    // Getting weights and max weight
    vector<double> weights;
    double maxWeight = std::numeric_limits<double>::min();
    for (int i = 0; i < num_particles; i++) {
        weights.push_back(particles[i].weight);
        if (particles[i].weight > maxWeight) {
            maxWeight = particles[i].weight;
        }
    }
    std::uniform_real_distribution<float> dist_float(0.0, maxWeight);
    std::uniform_real_distribution<float> dist_int(0.0, num_particles - 1);
    int index = dist_int(gen);
    double beta = 0.0;
    vector<Particle> resampledParticles;
    for (int i = 0; i < num_particles; i++) {
        beta += dist_float(gen) * 2.0;
        while (beta > weights[index]) {
            beta -= weights[index];
            index = (index + 1) % num_particles;
        }
        resampledParticles.push_back(particles[index]);
    }
    particles = resampledParticles;  

}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}
