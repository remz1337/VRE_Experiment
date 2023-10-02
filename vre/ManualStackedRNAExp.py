#Example command line arguments

import mysql.connector
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from argparse import ArgumentParser
from enum import Enum
import numpy as np
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms
from matplotlib import cm
from matplotlib.ticker import LinearLocator
from scipy.stats import chi2
import re
import glob
import os
import imageio
import math
import configparser
from collections import defaultdict
import graphviz
import traceback
from shutil import copyfile

from sqlalchemy import create_engine

from numpy import linalg as LA

#from scipy.ndimage.filters import maximum_filter
#from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
from scipy.ndimage import maximum_filter
from scipy.ndimage import generate_binary_structure, binary_erosion

from scipy import stats

import matplotlib.collections as mcoll
import matplotlib.path as mpath

from matplotlib.offsetbox import AnchoredText

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

#Uncomment this to fix matplotlib failling to allocate bitmap due to out of bound integer (happens on Windows with Visual studio)
#matplotlib.use('agg')

from matplotlib import rc
rc('font',**{'family':'serif','serif':['DejaVu Sans'],'size':15})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

#source:https://jonchar.net/notebooks/matplotlib-styling/
# Set the font used for MathJax - more on this later
rc('mathtext',**{'default':'regular'})


#Use pgf for final latex thesis, change to png for debugging. Handled by the --latex command line argument
#PLOT_EXTENSION=".pgf"

#Simply generate both
#PDF_EXTENSION=".pdf"
#PNG_EXTENSION=".png"

PLOT_FORMAT="pdf"
PLOT_EXTENSION="."+PLOT_FORMAT

BR_PLOT_Y_MAX=100

class Matter(object):
    pass

class ExperimentType(Enum):
    VRE = 1
    BE = 2
    BR = 3
    UNKNOWN = 4

def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument("-c", "--credentials", dest="credentials", help="path to MySQL config file containing the credentials.", metavar="FILE")
    parser.add_argument("-e", "--experiment", dest="experiment", help="Experiment ID to analyze and build the plots.")
    parser.add_argument("-f", "--folder", dest="folder", help="Root folder of the VRE experimentation. Used to save the plots along with the walks.")
    #parser.add_argument("-l", "--latex", dest="latex", help="Use this flag to save plots in LaTeX format, using PGF.", action='store_true')

    args = parser.parse_args()
    return args

def ValidateExperimentType(mysql_config,experiment_id):
    try:
        statement = "SELECT e.ID, e.type, e.track_variance, p.name, p.bnk_gates, p.bnk_inputs, a.generations, e.bnk_source_exp_id FROM experiment as e inner join problem as p on p.ID=e.fk_problem_id inner join algorithm as a on a.ID=e.fk_algorithm_id where e.ID = %(exp_id)s"
        parameters={"exp_id":experiment_id}

        engine_str="mysql+mysqlconnector://"+mysql_config["user"]+":"+mysql_config["password"]+"@"+mysql_config["host"]+":"+mysql_config["port"]+"/"+mysql_config["database"]
        engine = create_engine(engine_str, echo=False)
        df_record=pd.read_sql_query(statement, engine, params=parameters,)
    except Exception as e:
        print("Failed to read MySQL table")
        print("Caught exception:" + str(e))
        print(traceback.format_exc())

    return df_record

def RetrieveSamples(mysql_config,experiment_id):
    try:
        statement = "select ID from sample where fk_experiment_ID = %(exp_id)s"
        parameters={"exp_id":experiment_id}
        
        engine_str="mysql+mysqlconnector://"+mysql_config["user"]+":"+mysql_config["password"]+"@"+mysql_config["host"]+":"+mysql_config["port"]+"/"+mysql_config["database"]
        engine = create_engine(engine_str, echo=False)
        df_record=pd.read_sql_query(statement, engine, params=parameters,)
    except Exception as e:
        print("Failed to read MySQL table")
        print("Caught exception:" + str(e))
        print(traceback.format_exc())
    return df_record

def RetrieveWalks(mysql_config,sample_id):
    try:
        statement = "select ID, initial_population_size from walk where fk_sample_id = %(sample_id)s"
        parameters={"sample_id":sample_id}
        
        engine_str="mysql+mysqlconnector://"+mysql_config["user"]+":"+mysql_config["password"]+"@"+mysql_config["host"]+":"+mysql_config["port"]+"/"+mysql_config["database"]
        engine = create_engine(engine_str, echo=False)
        df_record=pd.read_sql_query(statement, engine, params=parameters,)
    except Exception as e:
        print("Failed to read MySQL table")
        print("Caught exception:" + str(e))
        print(traceback.format_exc())
    return df_record

def RetrieveVREWalkData(mysql_config,walk_id):
    try:
        statement = "select neighborhood_distance, local_genotype_robustness, local_genotype_evolvability, local_phenotype_robustness, local_phenotype_evolvability, total_neighborhood_size, total_neutral_neighbors, total_unique_phenotypes from walk_vre_stats where fk_walk_id = %(walk_id)s order by neighborhood_distance"
        parameters={"walk_id":walk_id}
        
        engine_str="mysql+mysqlconnector://"+mysql_config["user"]+":"+mysql_config["password"]+"@"+mysql_config["host"]+":"+mysql_config["port"]+"/"+mysql_config["database"]
        engine = create_engine(engine_str, echo=False)
        df_record=pd.read_sql_query(statement, engine, params=parameters,)
    except Exception as e:
        print("Failed to read MySQL table")
        print("Caught exception:" + str(e))
        print(traceback.format_exc())
    return df_record

def RetrieveVRESampleData(mysql_config,sample_id):
    try:
        statement = 'select neighborhood_distance,' \
            ' genotype_robustness_mean, genotype_evolvability_mean, phenotype_robustness_mean, phenotype_evolvability_mean,' \
            ' genotype_robustness_variance, genotype_evolvability_variance, phenotype_robustness_variance, phenotype_evolvability_variance' \
            ' from sample_vre_stats where fk_sample_id = %(sample_id)s order by neighborhood_distance'
        parameters={"sample_id":sample_id}
        
        engine_str="mysql+mysqlconnector://"+mysql_config["user"]+":"+mysql_config["password"]+"@"+mysql_config["host"]+":"+mysql_config["port"]+"/"+mysql_config["database"]
        engine = create_engine(engine_str, echo=False)
        df_record=pd.read_sql_query(statement, engine, params=parameters,)
    except Exception as e:
        print("Failed to read MySQL table")
        print("Caught exception:" + str(e))
        print(traceback.format_exc())
    return df_record

def RetrieveBRWalkData(mysql_config,walk_id):
    try:
        statement = 'select ID, combinations_count, combinations_genotype_distance_mean, combinations_phenotype_distance_mean,' \
            ' combinations_covariance_0_0,combinations_covariance_0_1,combinations_covariance_1_0,combinations_covariance_1_1,' \
            ' combinations_eigen_value_0,combinations_eigen_value_1,' \
            ' combinations_eigen_vector_0_0,combinations_eigen_vector_0_1,combinations_eigen_vector_1_0,combinations_eigen_vector_1_1, ellipse_angle, r_squared' \
            ' from baseline_robustness where fk_walk_id = %(walk_id)s'
        parameters={"walk_id":walk_id}
        
        engine_str="mysql+mysqlconnector://"+mysql_config["user"]+":"+mysql_config["password"]+"@"+mysql_config["host"]+":"+mysql_config["port"]+"/"+mysql_config["database"]
        engine = create_engine(engine_str, echo=False)
        df_record=pd.read_sql_query(statement, engine, params=parameters,)
    except Exception as e:
        print("Failed to read MySQL table")
        print("Caught exception:" + str(e))
        print(traceback.format_exc())
    return df_record

def RetrieveBRSampleData(mysql_config,sample_id):
    try:
        statement = 'select ID, genotype_distance_mean, phenotype_distance_mean, eigen_value_0_mean, eigen_value_1_mean, ellipse_angle_mean, r_squared_mean,' \
            ' genotype_distance_variance, phenotype_distance_variance, eigen_value_0_variance, eigen_value_1_variance, ellipse_angle_variance, r_squared_variance' \
            ' from sample_baseline_robustness where fk_sample_id = %(sample_id)s'
        parameters={"sample_id":sample_id}
        
        engine_str="mysql+mysqlconnector://"+mysql_config["user"]+":"+mysql_config["password"]+"@"+mysql_config["host"]+":"+mysql_config["port"]+"/"+mysql_config["database"]
        engine = create_engine(engine_str, echo=False)
        df_record=pd.read_sql_query(statement, engine, params=parameters,)
    except Exception as e:
        print("Failed to read MySQL table")
        print("Caught exception:" + str(e))
        print(traceback.format_exc())
    return df_record

def RetrieveBEWalkData(mysql_config,walk_id):
    try:
        statement = "select ID, value as baseline_value from baseline_evolvability where fk_walk_id = %(walk_id)s"
        parameters={"walk_id":walk_id}
        
        engine_str="mysql+mysqlconnector://"+mysql_config["user"]+":"+mysql_config["password"]+"@"+mysql_config["host"]+":"+mysql_config["port"]+"/"+mysql_config["database"]
        engine = create_engine(engine_str, echo=False)
        df_record=pd.read_sql_query(statement, engine, params=parameters,)
    except Exception as e:
        print("Failed to read MySQL table")
        print("Caught exception:" + str(e))
        print(traceback.format_exc())
    return df_record

def RetrieveBESampleData(mysql_config,sample_id):
    try:
        statement = 'select ID, baseline_evolvability_mean, baseline_evolvability_variance' \
            ' from sample_baseline_evolvability where fk_sample_id = %(sample_id)s'
        parameters={"sample_id":sample_id}
        
        engine_str="mysql+mysqlconnector://"+mysql_config["user"]+":"+mysql_config["password"]+"@"+mysql_config["host"]+":"+mysql_config["port"]+"/"+mysql_config["database"]
        engine = create_engine(engine_str, echo=False)
        df_record=pd.read_sql_query(statement, engine, params=parameters,)
    except Exception as e:
        print("Failed to read MySQL table")
        print("Caught exception:" + str(e))
        print(traceback.format_exc())
    return df_record


def PlotFitnessStats(walk_folder,BNK_gates,BNK_inputs):
    fd_title="BNK("+str(BNK_gates)+","+str(BNK_inputs)+") Fitness Distribution"
    #destination_file=walk_folder+"/"+fd_title+PLOT_EXTENSION
    destination_file=walk_folder+"/"+fd_title

    pop_file = walk_folder+"/vre.csv"

    #Number of bins
    num_bins=5

    #load data
    pop_data=pd.read_csv(pop_file).drop_duplicates()

    df2=pop_data.groupby(['Fitness']).size().reset_index(name='counts')

    df3=df2.sort_values("counts", ascending=False).head(num_bins)
    df3=df3.round({'Fitness': 6})

    df4 = df2.tail(-num_bins)

    df3.reset_index(drop=True, inplace=True)

    #df3.at[num_bins,'Fitness']='Other'
    df3.loc[num_bins,'Fitness']='Other'

    
    #df3.at[num_bins,'counts']=df4['counts'].sum()
    df3.loc[num_bins,'counts']=df4['counts'].sum()

    ax = df3.plot.bar(x='Fitness', y='counts', rot=0)

    plt.xticks(rotation = 45) # Rotates X-Axis Ticks by 45-degrees

    #plt.title(fd_title)

    plt.tight_layout()

    plt.savefig(destination_file+PLOT_EXTENSION)
    plt.close()

    return


def BuildExperimentFolderStr(root_folder, experiment_id):
    folder=root_folder+"/exp_"+str(experiment_id)
    return folder

def BuildSampleFolderStr(experiment_folder, sample_id):
    folder=experiment_folder+"/sample_"+str(sample_id)
    return folder

def BuildWalkFolderStr(sample_folder, walk_id):
    folder=sample_folder+"/walk_"+str(walk_id)
    return folder

def BuildNeighborhoodFolderStr(walk_folder, neighborhood_step):
    folder=walk_folder+"/neighborhood_"+str(neighborhood_step)
    return folder


def PlotFitness(walk_folder):
    pop_file = walk_folder+"/vre.csv"

    #load data
    pop_data=pd.read_csv(pop_file)

    title_bnk="Fitness evolution of a single walk"
    destination_file=walk_folder+"/Fitness"
    ax_bnk=plt.gca()

    #ax_bnk.set_title(title_bnk)
    ax_bnk.set_xlabel("Generation")
    ax_bnk.set_ylabel("Fitness")

    ax_bnk.plot(pop_data["Fitness"])

    plt.tight_layout()

    #plt.show()
    plt.savefig(destination_file+PLOT_EXTENSION)
    plt.close()

    return pop_data["Fitness"]


def PlotFitnessSample(df_fitness_data, sample_folder, max_generations):
    #destination_file=sample_folder+"/BE_plot"+PLOT_EXTENSION
    destination_file=sample_folder+"/Fitness"

    last_fitness=dict()

    best_fitness=[]
    mean_fitness=[]
    mean_fitness_var=[]

    for generation in range(max_generations):
        #print("generation "+str(generation))
        generation_fitness=[]
        for walk_id in df_fitness_data:
            fitness_val=-1 #init
            if(generation<len(df_fitness_data[walk_id])):
                fitness_val=df_fitness_data[walk_id][generation]
                last_fitness[walk_id]=fitness_val
            else:
                #terminating early means solution was found (ideal fitness)
                #we keep it to ensure it gets counted in the mean fitness across walks
                fitness_val=last_fitness[walk_id]
            #print(fitness_val)
            generation_fitness.append(fitness_val)

        best_fitness.append(max(generation_fitness))
        mean_fitness.append(np.mean(generation_fitness))
        mean_fitness_var.append(np.var(generation_fitness))


    x_count=range(1,max_generations+1)

    ax_be = plt.gca()

    #ax_be.plot(mean_fitness,label="mean")
    ax_be.errorbar(x_count, mean_fitness, yerr=mean_fitness_var, label='mean',  ecolor='silver')
    #ax_be.plot(best_fitness,label="best")
    #ax_be.legend()

    # Layout and titles
    #ax_be.set_title("Fitness evolution of a sample")
    ax_be.set_xlabel("Generation")
    ax_be.set_ylabel("Fitness")

    ##Add final Eb value on the graph
    ##props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
    #textstr="baseline evolvability = {0:.6f}".format(final_be)
    #text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5)
    ##at = AnchoredText("Figure 1a", loc='upper left', prop=dict(size=8), frameon=True)
    #text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    #ax_be.add_artist(text_box)

    #plt.setp(text_box.patch, facecolor='white', alpha=0.5)
    #ax_be.add_artist(text_box)

    #plt.show()
    plt.savefig(destination_file+PLOT_EXTENSION)
    plt.close()
    #https://timodenk.com/blog/exporting-matplotlib-plots-to-latex/

    #prepare dataframe for experiment plot

    df_result=pd.DataFrame()
    df_result["mean"]=mean_fitness
    df_result["best"]=best_fitness
    df_result["var"]=mean_fitness_var

    return df_result








def PlotStackedRNA(df_sample_30, df_sample_50, df_sample_100, base_folder):
    #destination_file=sample_folder+"/BE_plot"+PLOT_EXTENSION
    destination_file=base_folder+"/StackedRNA"

    fig=plt.figure()
    ax_rna=fig.gca()

    #ax_rna.errorbar(x_count, mean_fitness, yerr=mean_fitness_var, label='mean',  ecolor='silver')

    sample_id_30=next(iter(df_sample_30))
    sample_id_50=next(iter(df_sample_50))
    sample_id_100=next(iter(df_sample_100))

    #Assume only 1 sample!
    ax_rna.plot(df_sample_30[sample_id_30]["mean"],"--",label="L=30")
    ax_rna.plot(df_sample_50[sample_id_50]["mean"],"-.",label="L=50")
    ax_rna.plot(df_sample_100[sample_id_100]["mean"],":",label="L=100")
    # Layout and titles
    #ax_be.set_title("Fitness evolution of a sample")
    ax_rna.set_xlabel("generation")
    ax_rna.set_ylabel("mean fitness")

    ax_rna.legend()

    #plt.show()
    plt.savefig(destination_file+PLOT_EXTENSION)
    plt.close()

    return












def PlotFitnessExperiment(df_fitness_data, experiment_folder, max_generations, final_be=0):
    #destination_file=sample_folder+"/BE_plot"+PLOT_EXTENSION
    destination_file=experiment_folder+"/Fitness"

    best_fitness=[]
    mean_fitness=[]
    mean_fitness_var=[]

    for generation in range(max_generations):
        generation_mean_fitness=[]
        generation_best_fitness=[]
        for sample_id in df_fitness_data:
            #print(df_fitness_data[sample_id])
            #df_fitness_data[sample_id]["mean"]
            generation_mean_fitness.append(df_fitness_data[sample_id]["mean"].iloc[generation])
            generation_best_fitness.append(df_fitness_data[sample_id]["best"].iloc[generation])

        best_fitness.append(max(generation_best_fitness))
        mean_fitness.append(np.mean(generation_mean_fitness))
        mean_fitness_var.append(np.var(generation_mean_fitness))

    ax_be = plt.gca()

    x_count=range(1,max_generations+1)

    #ax_be.plot(mean_fitness,label="mean")
    ax_be.errorbar(x_count, mean_fitness, yerr=mean_fitness_var, label='mean',  ecolor='silver')
    #ax_be.plot(best_fitness,label="best")
    #ax_be.legend()

    # Layout and titles
    #ax_be.set_title("Fitness evolution of an experiment")
    ax_be.set_xlabel("Generation")
    ax_be.set_ylabel("Fitness")

    #Add final Eb value on the graph
    #props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
    textstr="baseline evolvability = {0:.6f}".format(final_be)
    text_box = AnchoredText(textstr, frameon=True, loc=4, pad=0.5)
    text_box.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
    #plt.setp(text_box.patch, facecolor='white', alpha=0.5)
    ax_be.add_artist(text_box)

    #plt.show()
    plt.savefig(destination_file+PLOT_EXTENSION)
    plt.close()
    #https://timodenk.com/blog/exporting-matplotlib-plots-to-latex/

    return


def main(args):
    #Maximum number of plots to generate (we don't need them only, only a few to show the results
    MAX_WALK_PLOTS=50#100
    MAX_SAMPLE_PLOTS=5#25

    global PLOT_EXTENSION
    credentials=args.credentials
    experiment_id=args.experiment
    root_folder=args.folder
    #use_pgf=args.latex

    exps = experiment_id.split(",")

    exp_id_30=int(exps[0])
    exp_id_50=int(exps[1])
    exp_id_100=int(exps[2])

    #if(use_pgf):
    #    PLOT_EXTENSION=".pgf"
    #else:
    #    PLOT_EXTENSION=".png"

    #Parse mysql credentials
    config = configparser.ConfigParser()
    config.read(credentials)

    mysql_config = {
      'user': config['client']['user'],
      'password': config['client']['password'],
      'host': config['client']['host'],
      'database': config['client']['database'],
      'port': config['client']['port']
    }

    #sample_measures=dict()
    #walk_measures=dict()
    sample_fitness_30=dict()
    sample_fitness_50=dict()
    sample_fitness_100=dict()
    walk_fitness_30=dict()
    walk_fitness_50=dict()
    walk_fitness_100=dict()
    #experiment_type=ExperimentType.UNKNOWN
    #track_variance=False

    experiment_30_folder=BuildExperimentFolderStr(root_folder,exp_id_30)
    experiment_50_folder=BuildExperimentFolderStr(root_folder,exp_id_50)
    experiment_100_folder=BuildExperimentFolderStr(root_folder,exp_id_100)

    experiment_details=ValidateExperimentType(mysql_config,exp_id_30)#Should be all the same anyway
    #exp_type_str=experiment_details["type"].iloc[0]
    #exp_track_variance_int=experiment_details["track_variance"].iloc[0]
    #problem_name_str=experiment_details["name"].iloc[0]
    #BNK_gates=experiment_details["bnk_gates"].iloc[0]
    #BNK_inputs=experiment_details["bnk_inputs"].iloc[0]
    max_generations=experiment_details["generations"].iloc[0]
    #bnk_source_exp_id_str=experiment_details["bnk_source_exp_id"].iloc[0]
    #bnk_source_exp_id=0
    #if bnk_source_exp_id_str is not None:
        #bnk_source_exp_id=int(bnk_source_exp_id_str)


    #generate_BNK_STD=False
    #Check if BNK
    #match = re.search('bnk', problem_name_str, re.IGNORECASE)
    #if match:
      #generate_BNK_STD=True

    #if exp_type_str == 'vre':
    #    experiment_type=ExperimentType.VRE
    #elif exp_type_str == 'br':
    #    experiment_type=ExperimentType.BR
    #elif exp_type_str == 'be':
    #    experiment_type=ExperimentType.BE

    #if exp_track_variance_int == 1:
    #    track_variance=True

    df_samples_30=RetrieveSamples(mysql_config,exp_id_30)
    df_samples_50=RetrieveSamples(mysql_config,exp_id_50)
    df_samples_100=RetrieveSamples(mysql_config,exp_id_100)

    #Retrieve all sample data
    #sample_cnt=0
    #for (sample_id_arr) in zip(df_samples_30["ID"]):
    for sample_id in df_samples_30["ID"]:
        #walk_cnt=0
        #sample_id=sample_id_arr[0]
        sample_folder=BuildSampleFolderStr(experiment_30_folder, sample_id)
        #if experiment_type == ExperimentType.VRE:
        #    sample_measures[sample_id]=RetrieveVRESampleData(mysql_config,sample_id)
        #elif experiment_type == ExperimentType.BR:
        #    sample_measures[sample_id]=RetrieveBRSampleData(mysql_config,sample_id)
        #elif experiment_type == ExperimentType.BE:
        #    sample_measures[sample_id]=RetrieveBESampleData(mysql_config,sample_id)

        df_walks=RetrieveWalks(mysql_config,sample_id)
        #for (walk_id_arr) in zip(df_walks["ID"]):
        for walk_id in df_walks["ID"]:
            #walk_id=walk_id_arr[0]
            walk_folder=sample_folder+"/walk_"+str(walk_id)

            walk_fitness_30[walk_id]=PlotFitness(walk_folder)                   

            #walk_cnt += 1
           
        #Build the plots for each sample
        sample_fitness_30[sample_id]=PlotFitnessSample(walk_fitness_30, sample_folder, max_generations)

        #sample_cnt += 1

    for sample_id in df_samples_50["ID"]:
        sample_folder=BuildSampleFolderStr(experiment_50_folder, sample_id)

        df_walks=RetrieveWalks(mysql_config,sample_id)
        for walk_id in df_walks["ID"]:
            walk_folder=sample_folder+"/walk_"+str(walk_id)
            walk_fitness_50[walk_id]=PlotFitness(walk_folder)                   
           
        #Build the plots for each sample
        sample_fitness_50[sample_id]=PlotFitnessSample(walk_fitness_50, sample_folder, max_generations)


    for sample_id in df_samples_100["ID"]:
        sample_folder=BuildSampleFolderStr(experiment_100_folder, sample_id)

        df_walks=RetrieveWalks(mysql_config,sample_id)
        for walk_id in df_walks["ID"]:
            walk_folder=sample_folder+"/walk_"+str(walk_id)
            walk_fitness_100[walk_id]=PlotFitness(walk_folder)                   
           
        #Build the plots for each sample
        sample_fitness_100[sample_id]=PlotFitnessSample(walk_fitness_100, sample_folder, max_generations)

    PlotStackedRNA(sample_fitness_30,sample_fitness_50,sample_fitness_100, experiment_30_folder)

    
    return

if __name__ == '__main__':
    args = parse_arguments()
    if (args.credentials != None):
        try:
            main(args)
        except Exception as e:
            print("Caught exception:" + str(e))
            print(traceback.format_exc())
    else:
        print("Missing arguments!")
