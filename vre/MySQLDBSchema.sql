CREATE TABLE `algorithm` (
  ID int NOT NULL AUTO_INCREMENT,
  name varchar(255) DEFAULT NULL,
  population_size int DEFAULT NULL,
  generations int DEFAULT NULL,
  other_parameters text DEFAULT NULL,
  PRIMARY KEY (ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE `problem` (
  ID int NOT NULL AUTO_INCREMENT,
  name varchar(255) DEFAULT NULL,
  genotype_length int DEFAULT NULL,
  rna_alphabet_size int DEFAULT NULL,
  bnk_gates int DEFAULT NULL,
  bnk_inputs int DEFAULT NULL,
  other_parameters text DEFAULT NULL,
  PRIMARY KEY (ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE `experiment` (
  ID int NOT NULL AUTO_INCREMENT,
  type varchar(255) DEFAULT NULL,
  comments varchar(255) DEFAULT NULL,
  max_samples int DEFAULT NULL,
  max_walks int DEFAULT NULL,
  target_phenotype text DEFAULT NULL,
  phenotype_distance_type VARCHAR(3) DEFAULT NULL;
  track_variance BOOLEAN NOT NULL DEFAULT FALSE,
  folder varchar(255) DEFAULT NULL,
  reruns int NOT NULL DEFAULT 0,
  neighborhoods int NOT NULL DEFAULT 1,
  bnk_source_exp_id int NULL DEFAULT NULL,
  fk_problem_id int DEFAULT NULL,
  fk_algorithm_id int DEFAULT NULL,
  PRIMARY KEY (ID),
  CONSTRAINT uc_experiment UNIQUE (fk_problem_id, fk_algorithm_id, max_samples, max_walks, reruns),
  FOREIGN KEY (fk_problem_id) REFERENCES problem(ID),
  FOREIGN KEY (fk_algorithm_id) REFERENCES algorithm(ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE `sample` (
  ID int NOT NULL AUTO_INCREMENT,
  random_id varchar(35) DEFAULT NULL,
  time_stamp timestamp DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  folder text DEFAULT NULL,
  already_analyzed BOOLEAN NOT NULL DEFAULT FALSE,
  fk_experiment_ID int DEFAULT NULL,
  PRIMARY KEY (ID),
  CONSTRAINT uc_sample UNIQUE (random_id),
  FOREIGN KEY (fk_experiment_ID) REFERENCES experiment(ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE `walk` (
  ID int NOT NULL AUTO_INCREMENT,
  random_id varchar(35) DEFAULT NULL,
  time_stamp timestamp DEFAULT CURRENT_TIMESTAMP ON UPDATE CURRENT_TIMESTAMP,
  initial_population_size int DEFAULT NULL,
  fk_sample_id int DEFAULT NULL,
  PRIMARY KEY (ID),
  CONSTRAINT uc_walk UNIQUE (random_id),
  FOREIGN KEY (fk_sample_id) REFERENCES sample(ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE `vre_neighborhood_job` (
  ID int NOT NULL AUTO_INCREMENT,
  distance int DEFAULT NULL,
  folder text DEFAULT NULL,
  java_args text,
  population_file text DEFAULT NULL,
  generated_neighbors int DEFAULT NULL,
  started BOOLEAN NOT NULL DEFAULT FALSE,
  started_time TIMESTAMP NULL DEFAULT NULL,
  slurm_job INT NULL DEFAULT NULL,
  completed BOOLEAN NOT NULL DEFAULT FALSE,
  completed_time TIMESTAMP NULL DEFAULT NULL,
  fk_walk_id int DEFAULT NULL,
  PRIMARY KEY (ID),
  FOREIGN KEY (fk_walk_id) REFERENCES walk(ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE `baseline_job` (
  ID int NOT NULL AUTO_INCREMENT,
  folder text DEFAULT NULL,
  java_args text,
  target_phenotype text DEFAULT NULL,
  population_file text DEFAULT NULL,
  started BOOLEAN NOT NULL DEFAULT FALSE,
  started_time TIMESTAMP NULL DEFAULT NULL,
  slurm_job INT NULL DEFAULT NULL,
  completed BOOLEAN NOT NULL DEFAULT FALSE,
  completed_time TIMESTAMP NULL DEFAULT NULL,
  fk_walk_id int DEFAULT NULL,
  PRIMARY KEY (ID),
  FOREIGN KEY (fk_walk_id) REFERENCES walk(ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE `sample_baseline_evolvability` (
  ID int NOT NULL AUTO_INCREMENT,
  baseline_evolvability_mean double DEFAULT NULL,
  baseline_evolvability_variance double DEFAULT NULL,
  fk_sample_ID int DEFAULT NULL,
  PRIMARY KEY (ID),
  CONSTRAINT uc_sample UNIQUE (fk_sample_ID),
  FOREIGN KEY (fk_sample_ID) REFERENCES sample(ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE `sample_baseline_robustness` (
  ID int NOT NULL AUTO_INCREMENT,
  genotype_distance_mean double DEFAULT NULL,
  phenotype_distance_mean double DEFAULT NULL,
  eigen_value_0_mean double DEFAULT NULL,
  eigen_value_1_mean double DEFAULT NULL,
  genotype_distance_variance double DEFAULT NULL,
  phenotype_distance_variance double DEFAULT NULL,
  eigen_value_0_variance double DEFAULT NULL,
  eigen_value_1_variance double DEFAULT NULL,
  ellipse_angle_mean double DEFAULT NULL,
  ellipse_angle_variance double DEFAULT NULL,
  r_squared_mean double DEFAULT NULL,
  r_squared_variance double DEFAULT NULL,
  aspect_ratio_mean double DEFAULT NULL,
  aspect_ratio_variance double DEFAULT NULL,
  orientation_mean double DEFAULT NULL,
  orientation_variance double DEFAULT NULL,
  sample_orientation double DEFAULT NULL,
  sample_aspect_ratio double DEFAULT NULL,
  fk_sample_ID int DEFAULT NULL,
  PRIMARY KEY (ID),
  CONSTRAINT uc_sample UNIQUE (fk_sample_ID),
  FOREIGN KEY (fk_sample_ID) REFERENCES sample(ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE `walk_vre_stats` (
  ID int NOT NULL AUTO_INCREMENT,
  neighborhood_distance int DEFAULT NULL,
  initial_phenotypes int DEFAULT NULL,
  total_neighborhood_size int DEFAULT NULL,
  total_neutral_neighbors int DEFAULT NULL,
  total_unique_phenotypes int DEFAULT NULL,
  local_genotype_robustness double DEFAULT NULL,
  local_genotype_evolvability double DEFAULT NULL,
  local_phenotype_robustness double DEFAULT NULL,
  local_phenotype_evolvability double DEFAULT NULL,
  fk_walk_id int DEFAULT NULL,
  PRIMARY KEY (ID),
  CONSTRAINT uc_walk UNIQUE (fk_walk_id, neighborhood_distance),
  FOREIGN KEY (fk_walk_id) REFERENCES walk(ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE `sample_vre_stats` (
  ID int NOT NULL AUTO_INCREMENT,
  neighborhood_distance int DEFAULT NULL,
  genotype_robustness_mean double DEFAULT NULL,
  genotype_evolvability_mean double DEFAULT NULL,
  phenotype_robustness_mean double DEFAULT NULL,
  phenotype_evolvability_mean double DEFAULT NULL,
  genotype_robustness_variance double DEFAULT NULL,
  genotype_evolvability_variance double DEFAULT NULL,
  phenotype_robustness_variance double DEFAULT NULL,
  phenotype_evolvability_variance double DEFAULT NULL,
  fk_sample_id int DEFAULT NULL,
  PRIMARY KEY (ID),
  CONSTRAINT uc_sample UNIQUE (fk_sample_id, neighborhood_distance),
  FOREIGN KEY (fk_sample_id) REFERENCES sample(ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE `baseline_evolvability` (
  ID int NOT NULL AUTO_INCREMENT,
  value double DEFAULT NULL,
  last_individual int DEFAULT NULL,
  final_fitness double DEFAULT NULL,
  fk_walk_id int DEFAULT NULL,
  PRIMARY KEY (ID),
  CONSTRAINT uc_walk UNIQUE (fk_walk_id),
  FOREIGN KEY (fk_walk_id) REFERENCES walk(ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE `baseline_robustness` (
  ID int NOT NULL AUTO_INCREMENT,
  combinations_count int DEFAULT NULL,
  combinations_genotype_distance_mean double DEFAULT NULL,
  combinations_phenotype_distance_mean double DEFAULT NULL,
  combinations_covariance_0_0 double DEFAULT NULL,
  combinations_covariance_0_1 double DEFAULT NULL,
  combinations_covariance_1_0 double DEFAULT NULL,
  combinations_covariance_1_1 double DEFAULT NULL,
  combinations_eigen_value_0 double DEFAULT NULL,
  combinations_eigen_value_1 double DEFAULT NULL,
  combinations_eigen_vector_0_0 double DEFAULT NULL,
  combinations_eigen_vector_0_1 double DEFAULT NULL,
  combinations_eigen_vector_1_0 double DEFAULT NULL,
  combinations_eigen_vector_1_1 double DEFAULT NULL,
  ellipse_angle double DEFAULT NULL,
  r_squared double DEFAULT NULL,
  aspect_ratio double DEFAULT NULL,
  orientation double DEFAULT NULL,
  fk_walk_id int DEFAULT NULL,
  PRIMARY KEY (ID),
  CONSTRAINT uc_walk UNIQUE (fk_walk_id),
  FOREIGN KEY (fk_walk_id) REFERENCES walk(ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;

CREATE TABLE `available_phenotypes` (
  ID int NOT NULL AUTO_INCREMENT,
  phenotype text DEFAULT NULL,
  length int DEFAULT NULL,
  used int NOT NULL DEFAULT 0,
  problem varchar(5),
  algorithm varchar(5),
  PRIMARY KEY (ID)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4;