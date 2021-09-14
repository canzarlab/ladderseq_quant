
#ifndef KALLISTO_EMALGORITHM_H
#define KALLISTO_EMALGORITHM_H

#include "common.h"
#include "KmerIndex.h"
#include "MinCollector.h"
#include "weights.h"

#include <algorithm>
#include <numeric>
#include <iostream>
#include <limits>
#include <vector>
#include <fstream>
#include <map>


// smallest weight we expect is ~10^-4
// on most machines, TOLERANCE should be 2.22045e-15
//const double TOLERANCE = std::numeric_limits<double>::epsilon() * 10;
//const double TOLERANCE = 1e-100;
const double TOLERANCE = std::numeric_limits<double>::denorm_min();

struct EMAlgorithm {
  // ecmap is the ecmap from KmerIndex
  // counts is vector from collector, with indices corresponding to ec ids
  // target_names is the target_names_ from collector
  // TODO: initialize alpha a bit more intelligently
  EMAlgorithm(const std::vector<int>& counts,
              const KmerIndex& index,
              const MinCollector& tc,
              const std::vector<double>& all_means,
              const ProgramOptions& opt) :
    index_(index),
    tc_(tc),
    num_trans_(index.target_names_.size()),
    ecmap_(index.ecmap),
    counts_(counts),
    target_names_(index.target_names_),
    post_bias_(4096,1.0),
    alpha_(num_trans_, 1.0/num_trans_), // uniform distribution over targets
    alpha_band_(num_trans_, 1.0/num_trans_), // uniform distribution over targets
    rho_(num_trans_, 0.0),
    rho_set_(false),
    all_fl_means(all_means),
    opt(opt)
  {
    assert(all_fl_means.size() == index_.target_lens_.size());
    eff_lens_ = calc_eff_lens(index_.target_lens_, all_fl_means);
    weight_map_ = calc_weights (tc_.counts, ecmap_, eff_lens_);

    assert(target_names_.size() == eff_lens_.size());
  }


  // Ladder Seq Changes: constructor taking n the different count files from each band
  EMAlgorithm(const std::vector<int>& counts,
              std::vector<MinCollector> tc_collection_,
              const KmerIndex& index,
              const MinCollector& tc,
              const std::vector<double>& all_means,
              const ProgramOptions& opt) :
    index_(index),
    tc_(tc),
    tc_collection_(tc_collection_),
    num_trans_(index.target_names_.size()),
    ecmap_(index.ecmap),
    counts_(counts),
    target_names_(index.target_names_),
    post_bias_(4096,1.0),
    alpha_(num_trans_, 1.0/num_trans_), // uniform distribution over targets
    alpha_band_(num_trans_, 1.0/num_trans_), // uniform distribution over targets
    rho_(num_trans_, 0.0),
    rho_set_(false),
    all_fl_means(all_means),
    opt(opt)
  {
    assert(all_fl_means.size() == index_.target_lens_.size());
    eff_lens_ = calc_eff_lens(index_.target_lens_, all_fl_means);
    std::vector<std::vector<int>> counts_band;
    for(int band = 0; band<7; band++){
        counts_band.push_back(tc_collection_[band].counts);
    }
    weight_map_band_= calc_weights_band(counts_band, ecmap_, eff_lens_);
    assert(target_names_.size() == eff_lens_.size());
  }

  // Ladder Seq Changes: constructor taking n the different count files from each band along with the band probabilities
  EMAlgorithm(const std::vector<int>& counts,
              std::vector<MinCollector> tc_collection_,
              const KmerIndex& index,
              const MinCollector& tc,
              const std::vector<double>& all_means,
              const ProgramOptions& opt,
              std::vector<std::vector<double>> beta_tb_) :
    index_(index),
    tc_(tc),
    tc_collection_(tc_collection_),
    num_trans_(index.target_names_.size()),
    ecmap_(index.ecmap),
    counts_(counts),
    target_names_(index.target_names_),
    post_bias_(4096,1.0),
    alpha_(num_trans_, 1.0/num_trans_), // uniform distribution over targets
    alpha_band_(num_trans_, 1.0/num_trans_), // uniform distribution over targets
    rho_(num_trans_, 0.0),
    rho_set_(false),
    all_fl_means(all_means),
    opt(opt),
    beta_tb_(beta_tb_)
  {
    assert(all_fl_means.size() == index_.target_lens_.size());
    eff_lens_ = calc_eff_lens(index_.target_lens_, all_fl_means);
    std::vector<std::vector<int>> counts_band;
    for(int band = 0; band<7; band++){
        counts_band.push_back(tc_collection_[band].counts);
    }
    weight_map_band_= calc_weights_band(counts_band, ecmap_, eff_lens_);
    assert(target_names_.size() == eff_lens_.size());
  }

  ~EMAlgorithm() {}

  //Ladder Seq Changes : Determines the length bucket of a transcript
  int determineLengtBucket(int t_it){
    int transcriptLength = index_.target_lens_[t_it] + 200;

    //For low resolution betas
    if(beta_tb_.size() < 15){
        //std::cout<<"Using low resolution betas!"<<std::endl;
        if(transcriptLength>0 && transcriptLength<=1000){
                return 0;
            }else if(transcriptLength>1000 && transcriptLength<=1500){
                return 1;
            }else if(transcriptLength>1500 && transcriptLength<=2000){
                return 2;
            }else if(transcriptLength>2000 && transcriptLength<=3000){
                return 3;
            }else if(transcriptLength>3000 && transcriptLength<=4000){
                return 4;
            }else if(transcriptLength>4000 && transcriptLength<=6000){
                return 5;
            }else if(transcriptLength>6000){
                return 6;
            }
    }else{
    // High resolution betas -- new cuts
        if(transcriptLength>0 && transcriptLength<=910){
            return 0;
        }else if(transcriptLength>910 && transcriptLength<=1083){
            return 1;
        }else if(transcriptLength>1083 && transcriptLength<=1268){
            return 2;
        }else if(transcriptLength>1268 && transcriptLength<=1440){
            return 3;
        }else if(transcriptLength>1440 && transcriptLength<=1590){
            return 4;
        }else if(transcriptLength>1590 && transcriptLength<=1757){
            return 5;
        }else if(transcriptLength>1757 && transcriptLength<=1920){
            return 6;
        }else if(transcriptLength>1920 && transcriptLength<=2070){
            return 7;
        }else if(transcriptLength>2070 && transcriptLength<=2273){
            return 8;
        }else if(transcriptLength>2273 && transcriptLength<=2462){
            return 9;
        }else if(transcriptLength>2462 && transcriptLength<=2701){
            return 10;
        }else if(transcriptLength>2701 && transcriptLength<=2932){
            return 11;
        }else if(transcriptLength>2932 && transcriptLength<=3230){
            return 12;
        }else if(transcriptLength>3230 && transcriptLength<=3554){
            return 13;
        }else if(transcriptLength>3554 && transcriptLength<=3893){
            return 14;
        }else if(transcriptLength>3893 && transcriptLength<=4525){
            return 15;
        }else if(transcriptLength>4525 && transcriptLength<=5324){
            return 16;
        }else if(transcriptLength>5324 && transcriptLength<=6378){
            return 17;
        }else if(transcriptLength>6378){
            return 18;
        }

    }
  }

  double round(double val){
      if(val < 0) return ceil(val - 0.5);
      return floor(val + 0.5);
  }

  double sum(std::vector<double> givenVector){
      double sum = 0;
      for(auto a : givenVector){
        sum += a;
      }
      return sum;
  }
//Ladder seq changes end

  void run(size_t n_iter = 1000, size_t min_rounds=50, bool verbose = true, bool recomputeEffLen = true) {

    std::vector<double> next_alpha(alpha_.size(), 0.0);

    // Ladder Seq Changes: Making the next_alpha_tb vector //
    std::vector<std::vector<double>> next_alpha_tb(alpha_.size());
    for(int i= 0; i<next_alpha_tb.size(); i++){
        next_alpha_tb[i] = std::vector<double> (7,0.0);
    }
    // Ladder Seq Changes End

    //Reading in the band probabilities from the supplied file Only if the file name is provided
    if(opt.bandProbabilityFile.size()!=0){
        std::ifstream bandProbabilityFile(opt.bandProbabilityFile);
        std::string line;
        while (std::getline(bandProbabilityFile, line)){
            if(line.size()==0)continue;
            std::vector<double> beta_length;
            std::istringstream iss(line);
            double a = 0;
            while(iss>>a){
                beta_length.emplace_back(a);
            }
            beta_tb_.push_back(beta_length);
        }
        bandProbabilityFile.close();
    }
    // Ladder seq changes end //



//    assert(weight_map_.size() <= counts_.size());
    assert(weight_map_band_.size() <= counts_.size());

    double denom;
    std::vector<double> denomBand(7); //Ladder Seq Changes
    double denomBandTotal; //Ladder Seq Changes
    const double alpha_limit = 1e-7;
    const double alpha_change_limit = 1e-2;
    const double alpha_change = 1e-2;
    bool finalRound = false;

    if (verbose) {
      std::cerr << "[   em] quantifying the abundances ..."; std::cerr.flush();
    }

    // This vector would contain the individual counts from each band for each transcript
    // Format : transcript- Band1Count,Band2Count ...
    std::vector<std::vector<double>> individialBandCounts(alpha_.size());
    for(int i= 0; i<individialBandCounts.size(); i++){
        individialBandCounts[i] = std::vector<double> (7,0.0);
    }


//   std::cout<<"Hereeee111"<<std::endl;
    int i;
    for (i = 0; i < n_iter; ++i) {
      double totalNumberofReads = 0;

      if (recomputeEffLen && (i == min_rounds || i == min_rounds + 500)) {
        eff_lens_ = update_eff_lens(all_fl_means, tc_, index_, alpha_, eff_lens_, post_bias_, opt);
//        weight_map_ = calc_weights (tc_.counts, ecmap_, eff_lens_);

        //Ladder Seq Changes
        std::vector<std::vector<int>> counts_band;
        for(int band = 0; band<7; band++){
            counts_band.push_back(tc_collection_[band].counts);
        }
        weight_map_band_= calc_weights_band(counts_band, ecmap_, eff_lens_);
        //Ladder Seq Changes End
      }
//     std::cout<<"Hereeee222"<<std::endl;

      // Ladder Comment: loop through all the singular equivalence classes
      for (int ec = 0; ec < num_trans_; ec++) {
//        //LadderSeq Changes
          int sumOfCounts = 0;
        for(int band = 0; band<7; band++){
            (next_alpha_tb.at(ec)).at(band) = 0;
            (individialBandCounts.at(ec)).at(band) = tc_collection_[band].counts.at(ec);
            sumOfCounts += tc_collection_[band].counts.at(ec);
        }
        next_alpha[ec] = sumOfCounts;
        totalNumberofReads += sumOfCounts;

//        //LadderSeq Changes End
      }

//     std::cout<<"Hereeee333"<<std::endl;
      // Ladder Comment: loop through all the non singular equivalence classes
      for (int ec = num_trans_; ec < ecmap_.size();  ec++) {
        denom = 0.0;
        denomBandTotal = 0.0;
        std::fill(denomBand.begin(), denomBand.end(), 0.0); // Ladder Seq changes

        auto& v = ecmap_[ec]; //ecmap_.find(ec)->second;

        //LadderSeq Changes
        // Checking if there are counts in at least one band for a particular equivalence class
        bool flag = false;
        for(int band = 0; band<7; band++){
            int ec_index = tc_collection_[band].findEC(v);
            if(tc_collection_[band].counts.at(ec_index)!=0){
                flag = true;
            }
        }
        if(flag==false){
            continue;
        }
        //LadderSeq Changes End

//       std::cout<<"Hereeee444"<<std::endl;

        auto numEC = v.size();


        //Ladder Seq//
        //Calculating a different denominator for each band
        auto& wv_band = weight_map_band_[ec];

        for(int band = 0; band<7 ; band++){
            for (auto t_it = 0; t_it < numEC; ++t_it) {
                int lengthBucket = determineLengtBucket(v[t_it]);
              denomBand[band] += alpha_[v[t_it]] * wv_band[band][t_it] * beta_tb_[lengthBucket].at(band);
            }
            denomBandTotal += denomBand[band] ;
        }
        if (denomBandTotal < TOLERANCE) {
          continue;
        }

//       std::cout<<"Hereeee555"<<std::endl;
        // Compute the countNorm for each band
        std::vector<double> countNormBand(7,0.0);
        for(int band = 0; band<7; band++){
            if(denomBand[band]==0){
                countNormBand[band] = 0;
                continue;
            }
            int ec_index = tc_collection_[band].findEC(v);
            countNormBand[band] = (double)tc_collection_[band].counts.at(ec_index) / denomBand[band];
        }
//       std::cout<<"Hereeee666"<<std::endl;
//        sd::cout<<numEC<<std::endl;
        // compute the update step
        // First sum up all the different alphatb
        for(int band = 0; band<7; band++){
            for (auto t_it = 0; t_it < numEC; ++t_it) {
                if(wv_band[band][t_it]==0 || countNormBand[band]==0){
                    (next_alpha_tb.at(v[t_it])).at(band) += 0;
                    continue;
                }
                int lengthBucket = determineLengtBucket(v[t_it]);
                (next_alpha_tb.at(v[t_it])).at(band) +=  ((beta_tb_[lengthBucket]).at(band)) * wv_band[band][t_it] * alpha_[v[t_it]] * countNormBand[band];
            }
        }

        // Second sum up the different alphatb from each band
        for (auto t_it = 0; t_it < numEC; ++t_it) {
            for(int band = 0; band<7; band++){
                next_alpha[v[t_it]] += (next_alpha_tb.at(v[t_it])).at(band);
                (individialBandCounts.at(v[t_it])).at(band) += (next_alpha_tb.at(v[t_it])).at(band);
                (next_alpha_tb.at(v[t_it])).at(band) = 0;
            }
        }
//       std::cout<<"Hereeee777"<<std::endl;
        //Ladder Seq Ends//
      }


      // TODO: check for relative difference for convergence in EM

      bool stopEM = false; //!finalRound && (i >= min_rounds); // false initially
      //double maxChange = 0.0;
      int chcount = 0;


      for (int ec = 0; ec < num_trans_; ec++) {
        if (next_alpha[ec] > alpha_change_limit && (std::fabs(next_alpha[ec] - alpha_[ec]) / next_alpha[ec]) > alpha_change) {
          chcount++;
        }
        alpha_[ec] = next_alpha[ec];
        // clear all next_alpha values 0 for next iteration
        next_alpha[ec] = 0.0;
        //LadderSeq Changes: Clear all the next_alpha_tb values
        for(int band = 0; band<7; band++){
            (next_alpha_tb.at(ec)).at(band)= 0.0;
        }
      }



      //LadderSeq Changes End

      //std::cout << chcount << std::endl;
      if (chcount == 0 && i > min_rounds) {

        stopEM=true;
      }

      if (finalRound) {
        break;
      }

      // std::cout << maxChange << std::endl;
      if (stopEM) {
        finalRound = true;
        alpha_before_zeroes_.resize( alpha_.size() );
        for (int ec = 0; ec < num_trans_; ec++) {
          alpha_before_zeroes_[ec] = alpha_[ec];
          if (alpha_[ec] < alpha_limit/10.0) {
            alpha_[ec] = 0.0;
          }
        }
      }
    } // The Em algorithm ends here

    std::cout<<"Total alpha after EM: "<<this->sum(alpha_)<<std::endl;

     // Writing the fractional counts for each transcript
     std::ofstream transcriptCountFile;
     transcriptCountFile.open (opt.output+"/transcriptCountFile.txt");
     for (int ec = 0; ec < num_trans_; ec++){
    transcriptCountFile <<target_names_.at(ec)<<"\t";
        for (int band = 0; band<7; band++){
                transcriptCountFile<<individialBandCounts.at(ec).at(band)<<"\t";
            }
        transcriptCountFile<<"\n";
     }
     transcriptCountFile.close();

    // ran for the maximum number of iterations
    if (n_iter == i) {
      alpha_before_zeroes_.resize( alpha_.size() );
      for (int ec = 0; ec < num_trans_; ec++) {
        alpha_before_zeroes_[ec] = alpha_[ec];
      }
    }

    if (verbose) {
      std::cerr << " done" << std::endl;
      std::cerr << "[   em] the Expectation-Maximization algorithm ran for "
        << pretty_num(i) << " rounds";
      std::cerr << std::endl;
      std::cerr.flush();
    }
  }

  void compute_rho() {
    if (rho_set_) {
      // rho has already been set, let's clear it
      std::fill(rho_.begin(), rho_.end(), 0.0);
    }

    double total {0.0};
    for (auto i = 0; i < alpha_.size(); ++i) {
      if (eff_lens_[i] < TOLERANCE) {
        std::cerr << "Should actually never really get here... tid: "  << i <<
            std::endl;
        continue;
      }
      rho_[i] = alpha_[i] / eff_lens_[i];
      total += rho_[i];
    }

    for (auto& r : rho_) {
      r /= total;
    }

    rho_set_ = true;
  }

  // DEPRECATED:
  void write(const std::string& out_fname) const {
    std::ofstream out;
    out.open(out_fname, std::ios::out);

    if (!out.is_open()) {
      std::cerr << "Error opening '" << out_fname << "'" <<
          std::endl;
      exit(1);
    }

    out.precision(15);

    out <<
        "target_id" << "\t" <<
        "kallisto_id" << "\t" <<
        "rho" << "\t" <<
        "tpm" << "\t" <<
        "est_counts" <<
        std::endl;

    const double MILLION = 1e6;

    for (auto i = 0; i < rho_.size(); ++i) {
      out <<
          target_names_[i] << "\t" <<
          i << "\t" <<
          rho_[i] << "\t" <<
          rho_[i] * MILLION << "\t" <<
          alpha_[i] <<
          std::endl;
    }

    out.flush();
    out.close();
  }

  void set_start(const EMAlgorithm& em_start) {
    assert(em_start.alpha_before_zeroes_.size() == alpha_.size());
    double big = 1.0;
    double sum_counts = std::accumulate(counts_.begin(), counts_.end(), 0.0);
    double sum_big = 0.0;
    int count_big = 0;
    for (auto x : em_start.alpha_before_zeroes_) {
      if (x >= big) {
        sum_big += x;
        count_big++;
      }
    }
    int n = alpha_.size();
    for (auto i = 0; i < n; i++) {
      if (em_start.alpha_before_zeroes_[i] >= big) {
        alpha_[i] = em_start.alpha_before_zeroes_[i];
      } else {
        alpha_[i] = sum_counts/(n - count_big);
      }
    }

    //std::cout << sum_big << " " << count_big << " " << n << std::endl;

    std::copy(em_start.alpha_before_zeroes_.begin(), em_start.alpha_before_zeroes_.end(),
        alpha_.begin());
  }


  int num_trans_;
  const KmerIndex& index_;
  const MinCollector& tc_;
  const EcMap& ecmap_;
  const std::vector<int>& counts_;
  const std::vector<std::string>& target_names_;
  const std::vector<double>& all_fl_means;
  std::vector<double> eff_lens_;
  std::vector<double> post_bias_;
  WeightMap weight_map_;
  std::vector<double> alpha_;
  std::vector<double> alpha_before_zeroes_;
  std::vector<double> rho_;
  bool rho_set_;
  const ProgramOptions& opt;
  //Ladder Seq Changes//
//  std::vector<std::vector<double>> alpha_tb_;
  std::vector<std::vector<double>> beta_tb_;
  std::vector<MinCollector> tc_collection_;
  std::vector<WeightMap> weight_map_band_;
  std::vector<double> alpha_band_;

  //Ladder Seq Changes End//
};


#endif // KALLISTO_EMALGORITHM_H
