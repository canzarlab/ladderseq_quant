
#ifndef KALLISTO_MIGRATIONPROBABILITYESTIMATOR_H
#define KALLISTO_MIGRATIONPROBABILITYESTIMATOR_H

#include "common.h"
#include "KmerIndex.h"
#include "MinCollector.h"
#include "weights.h"



struct MigrationProbabilityEstimator {

  // Ladder Seq Changes: constructor taking n the different count files from each band
  MigrationProbabilityEstimator(std::vector<MinCollector> tc_collection_,
              const KmerIndex& index) :
    tc_collection_(tc_collection_),
    index_(index),
    num_trans_(index.target_names_.size()){
  }

  ~MigrationProbabilityEstimator() {}

  // Determines the length bucket of a transcript
  int determineLengtBucket(int t_it){
    int transcriptLength = index_.target_lens_[t_it] + 200;
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


  std::vector<std::vector<double>> estimate(bool verbose = true) {


    if (verbose) {
      std::cerr << "Estimating migration probabilities ... "; std::cerr.flush();
    }

    std::vector<std::vector<double>> migrationProbabilities(19);
    for(int i= 0; i<migrationProbabilities.size(); i++){
        migrationProbabilities[i] = std::vector<double> (7,0.0);
    }

    // Loop through all the singular equivalence classes
    for (int ec = 0; ec < num_trans_; ec++) {
      int sumOfCounts = 0;
      int lengthBucket = determineLengtBucket(ec);
      for(int band = 0; band<7; band++){
          sumOfCounts += tc_collection_[band].counts.at(ec);
      }
      if(sumOfCounts<51){
          continue;
      }
      for(int band = 0; band<7; band++){
          (migrationProbabilities.at(lengthBucket)).at(band) += tc_collection_[band].counts.at(ec);
      }
    }

    for (int i = 0; i < 19; i++) {
        int sum = 0;
        for(int band = 0; band<7; band++){
            sum += migrationProbabilities.at(i).at(band);
        }
        for(int band = 0; band<7; band++){
            migrationProbabilities.at(i).at(band) = migrationProbabilities.at(i).at(band)/sum;
//            std::cout<<migrationProbabilities.at(i).at(band)<<" ";
        }
//        std::cout<<" "<<std::endl;
    }


    return migrationProbabilities;
}



  std::vector<MinCollector> tc_collection_;
  int num_trans_;
  const KmerIndex& index_;
};


#endif // KALLISTO_EMALGORITHM_H
