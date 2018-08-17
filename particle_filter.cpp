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
#include <map>

#include "particle_filter.h"
#include "helper_functions.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	//cout << "init starts" << endl;

	num_particles = 100;
	weights.resize(num_particles);
	particles.resize(num_particles);

	default_random_engine genor;
	// normal distribution for sensor nosie
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

	 for (int i = 0; i < num_particles; ++i) {
	    particles[i].id = i;
	    particles[i].x = dist_x(genor);
	    particles[i].y = dist_y(genor);
	    particles[i].theta = dist_theta(genor);
	    particles[i].weight = 1.0;
	  }

	is_initialized = true;
	//cout << "init ends" << endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	//cout << "prediction starts" << endl;

	//implementation of equation in Lesson 14 to update x, y and yaw angle
	default_random_engine genor;
	normal_distribution<double>noiseX(0, std_pos[0]);
	normal_distribution<double>noiseY(0, std_pos[1]);
	normal_distribution<double>noiseTheta(0, std_pos[2]);

		for (int i=0; i<num_particles; ++i){
			if (abs(yaw_rate) > 0.0001){
				particles[i].x += (velocity/yaw_rate)*(sin(particles[i].theta + (yaw_rate*delta_t)) - sin(particles[i].theta)) ;
				particles[i].y += (velocity/yaw_rate)*(cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
				particles[i].theta += yaw_rate * delta_t;
			}else{
				particles[i].x += velocity*delta_t*cos(particles[i].theta);
				particles[i].y += velocity*delta_t*sin(particles[i].theta);
			}
			particles[i].x += noiseX(genor);
			particles[i].y += noiseY(genor);
			particles[i].theta += noiseTheta(genor);
		}

	//cout << "prediction ends" << endl;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

	//cout << "dataAssociation starts" << endl;

	//Loop over observation particles
	int predID;
	for (int i=0; i<observations.size();++i){
		LandmarkObs obs1 = observations[i];
		double minDist = INFINITY;
		int obsID = -1;
		for (int j=0;j<predicted.size();++j){
			LandmarkObs pred1 = predicted[j];
			double sumDist;
			sumDist = dist(obs1.x, obs1.y, pred1.x, pred1.y);
			if (sumDist < minDist){
				minDist = sumDist;
				obsID = i;
			}
		}
		observations[i].id = obsID;
		}
	//cout << "dataAssociation ends" << endl;
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
	//cout << "updateWeights starts" << endl;
	//construct multivariate Gaussian distribution:
	double term1 = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]);
	double termx = 2.0 * std_landmark[0] * std_landmark[0];
	double termy = 2.0 * std_landmark[1] * std_landmark[1];
	//cout << "updateWeights: term1 = " << term1 << endl;
	//cout << "updateWeights: termx = " << termx << endl;
	//cout << "updateWeights: termy = " << termy << endl;
	for (int i=0;i<num_particles;++i){
		//init
		double w = 1.0;
		for (int j=0;j<observations.size();++j){
			//transform observation coord to map coord
			double obsX, obsY;
			obsX =  observations[j].x*cos(particles[i].theta) - observations[j].y*sin(particles[i].theta) + particles[i].x;
			obsY =  observations[j].x*sin(particles[i].theta) + observations[j].y*cos(particles[i].theta) + particles[i].y;

			vector<Map::single_landmark_s> landmarks = map_landmarks.landmark_list;
			vector<double> landmarkObsDist (landmarks.size());
			for (int k=0;k<landmarks.size();++k){
				double landmarkParDist = dist(landmarks[k].x_f,landmarks[k].y_f,particles[i].x,particles[i].y);
				if (landmarkParDist<= sensor_range){
					landmarkObsDist[k] = dist(landmarks[k].x_f,landmarks[k].y_f,obsX,obsY);
				}else{
					landmarkObsDist[k] = numeric_limits<float>::max();
				}
				//cout << "updateWeights: loop over k = " << k << endl;
			}
			// find observation point
			int localPos =  distance(landmarkObsDist.begin(),min_element(landmarkObsDist.begin(),landmarkObsDist.end()));
			//cout << "updateWeights: localPos = " << localPos << endl;
			double muX = landmarks[localPos].x_f;
			double muY = landmarks[localPos].y_f;
			//cout << "updateWeights: muX = " << muX << endl;
			//cout << "updateWeights: muY = " << muY << endl;
			//cout << "updateWeights: obsX = " << obsX << endl;
			//cout << "updateWeights: obsY = " << obsY << endl;
			double xDiff = obsX - muX;
			double yDiff = obsY - muY;
			//cout << "updateWeights: xDiff = " << xDiff << endl;
			//cout << "updateWeights: yDiff = " << yDiff << endl;
			double tmp1 = (xDiff*xDiff)/termx + (yDiff*yDiff)/termy;
			//cout << "updateWeights: tmp1 = " << tmp1 << endl;
			double term4 = exp(-tmp1);
			//cout << "updateWeights: term4 = " << term4 << endl;
			w *= term1 * term4;
			//cout << "updateWeights: loop over j = " << j << endl;
		}
		//cout << "updateWeights: loop over i = " << i << endl;

		particles[i].weight = w;
		weights[i] = particles[i].weight;
		//cout << "updateWeights: w = " << w << endl;
		//cout << "updateWeights: weights = " << weights[i] << endl;
	}

	//cout << "updateWeights ends" << endl;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	//cout << "resample starts" << endl;

	default_random_engine genor;
	vector<Particle> newParticles(num_particles);

	for (int i=0;i<num_particles;++i){
		discrete_distribution<int>index(weights.begin(), weights.end());
		newParticles[i] = particles[index(genor)];
	}
	particles = newParticles;

	//cout << "resample ends" << endl;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

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
