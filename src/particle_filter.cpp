/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  num_particles = 100;
  particles.resize(num_particles);
  weights.resize(num_particles);

  for (int i = 0; i < num_particles; ++i)
  {
    particles[i].id = i;
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
    particles[i].weight = 1;
  }

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

  default_random_engine gen;

  // Create distribution around origin so we don't have to recalculate for every particle
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  //Non-zero yaw rate equations of motion with Gaussian noise
  if (fabs(yaw_rate) > 0.0001)
  {
    for (int i = 0; i < num_particles; ++i)
    {
      particles[i].x += (velocity / yaw_rate) * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)) + dist_x(gen);
      particles[i].y += (velocity / yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t)) + dist_y(gen);
      particles[i].theta += yaw_rate * delta_t + dist_theta(gen);
    }
  }
  // Zero yaw rate equations of motion with Gaussian noise
  else
  {
    for (int i = 0; i < num_particles; ++i)
    {
      particles[i].x += velocity * delta_t * cos(particles[i].theta) + dist_x(gen);
      particles[i].y += velocity * delta_t * sin(particles[i].theta) + dist_y(gen);
      particles[i].theta += dist_theta(gen);
    }
  }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

  double p_x, p_y, p_theta, obs_x, obs_y, t_obs_x, t_obs_y;
  double MGPD_norm = 1 / (2 * M_PI * std_landmark[0] * std_landmark[1]);

  for (int i = 0; i < num_particles; ++i)
  {
    p_x = particles[i].x;
    p_y = particles[i].y;
    p_theta = particles[i].theta;

    particles[i].weight = 1;  //  reset the weight for later calculation


    vector<LandmarkObs> observable_landmarks; // list of all observable landmarks

    //  find all landmarks within sensor range of current particle and append to list
    for (int j = 0; j < map_landmarks.landmark_list.size(); ++j)
    {
      LandmarkObs current_landmark = { map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f };

      // calculate distance from particle position to landmarks
      double landmark_dist = dist(p_x, p_y, current_landmark.x, current_landmark.y);

      if (landmark_dist <= sensor_range)
      {
        observable_landmarks.push_back(current_landmark);
      }
    }


    //  transform all observations from vehicle coordinates (local) to particle coordinates (global)
    for (int k = 0; k < observations.size(); ++k)
    {
        obs_x = observations[k].x;
        obs_y = observations[k].y;

        // transform to global coordinates
        t_obs_x = p_x + cos(p_theta) * obs_x - sin(p_theta) * obs_y;
        t_obs_y = p_y + sin(p_theta) * obs_x + cos(p_theta) * obs_y;

        //  find the landmark that has a minimum distance to the observation and assign the landmark ID to the observation ID
        double obs_dist = dist(t_obs_x, t_obs_y, observable_landmarks[0].x, observable_landmarks[0].y);
        double min_obs_x, min_obs_y;

        //  iterate through all observable landmarks
        for (int n = 0; n < observable_landmarks.size(); ++n)
        {
          double temp_obs_dist = dist(t_obs_x, t_obs_y, observable_landmarks[n].x, observable_landmarks[n].y);
          if (temp_obs_dist <= obs_dist)
          {
            obs_dist = temp_obs_dist;
            min_obs_x = observable_landmarks[n].x;
            min_obs_y = observable_landmarks[n].y;
          }
        }
    
        //  calculate the particle weight
        double x_exp = (min_obs_x - t_obs_x) * (min_obs_x - t_obs_x) / (2 * std_landmark[0] * std_landmark[0]);
        double y_exp = (min_obs_y - t_obs_y) * (min_obs_y - t_obs_y) / (2 * std_landmark[1] * std_landmark[1]);

        particles[i].weight *= MGPD_norm * exp(-(x_exp + y_exp));
    }
    weights[i] = particles[i].weight; //  update the weight vector
  }
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

  std::vector<Particle> resampled_particles;

  //  create a new distribution based on the weights vector
  std::discrete_distribution<int> distribution(weights.begin(),weights.end());

  //  create a new particle set based on the weight distribution
  for (int i = 0; i < num_particles; ++i)
  {
    resampled_particles.push_back(particles[distribution(gen)]);
  }

  particles = resampled_particles;  //  push the resampled particles to the particle list
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;

    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
