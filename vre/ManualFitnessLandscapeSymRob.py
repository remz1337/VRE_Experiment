#Example command line arguments

from curses import qiflush
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

import matplotlib.collections as mcoll
import matplotlib.path as mpath

from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

#Uncomment this to fix matplotlib failling to allocate bitmap due to out of bound integer (happens on Windows with Visual studio)
#matplotlib.use('agg')

from matplotlib import rc
rc('font',**{'family':'serif','serif':['DejaVu Sans']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

#source:https://jonchar.net/notebooks/matplotlib-styling/
# Set the font used for MathJax - more on this later
rc('mathtext',**{'default':'regular'})

#Use pgf for final latex thesis, change to png for debugging. Handled by the --latex command line argument
#PLOT_EXTENSION=".pgf"

#Simply generate both
PGF_EXTENSION=".pgf"
PNG_EXTENSION=".png"

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
    parser.add_argument("-l", "--latex", dest="latex", help="Use this flag to save plots in LaTeX format, using PGF.", action='store_true')

    args = parser.parse_args()
    return args

def ValidateExperimentType(mysql_config,experiment_id):
    try:
        statement = "SELECT e.ID, e.type, e.track_variance, p.name, p.bnk_gates, p.bnk_inputs, a.generations FROM experiment as e inner join problem as p on p.ID=e.fk_problem_id inner join algorithm as a on a.ID=e.fk_algorithm_id where e.ID = %(exp_id)s"
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
            ' combinations_eigen_vector_0_0,combinations_eigen_vector_0_1,combinations_eigen_vector_1_0,combinations_eigen_vector_1_1, ellipse_angle' \
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
        statement = 'select ID, genotype_distance_mean, phenotype_distance_mean, eigen_value_0_mean, eigen_value_1_mean, ellipse_angle_mean,' \
            ' genotype_distance_variance, phenotype_distance_variance, eigen_value_0_variance, eigen_value_1_variance, ellipse_angle_variance' \
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

def PlotBNKWalkSTD(pop_folder,BNK_gates,BNK_inputs):
    #https://stackoverflow.com/questions/47153044/python-transitions-how-to-map-all-possible-transitions-of-an-fsm-to-a-graph
    
    #First read pop file
    #For X individuals:
        #Parse genotype
        #Evaluate states
        #Build transition array for the transition library

    destination_folder=pop_folder+"/transition_diagrams"
    try:
        if not os.path.exists(destination_folder):
            os.mkdir(destination_folder)
    except OSError:
        print ("Creation of the directory %s failed" % destination_folder)

    pop_file = pop_folder+"/vre.csv"

    #load data
    pop_data=pd.read_csv(pop_file)

    output_len=pow(2, BNK_inputs)
    genes_per_gate = BNK_inputs+output_len

    TOTAL_STATES=pow(2,BNK_gates)
    format_rule='0'+ str(BNK_gates) +'b'

    #Maybe filter by fitness? or use the phenotype (-1) for the irregular ones? but need to limit the number...
    #for index, row in pop_data[pop_data['Phenotype'].str.contains("-1")].head(5).iterrows():
    #Can't know anymore if there are irregular cycles (and maybe we dont care since our new measure are able to handle them...)
    #so just plot the 5 first
    for index, row in pop_data.head(5).iterrows():
        genotype=row["Genotype"]
        fitness=row["Fitness"]
        phenotype=row["Phenotype"]

        individual_line_num=pop_data.index[pop_data['Genotype'] == genotype].tolist()[0]

        #Split by white space
        genes=genotype.split()

        state_diagram=list()
        states=list()

        for state in range(TOTAL_STATES):
            machine_state_str=format(state, format_rule)
            states.append(machine_state_str)
            next_state_str=""
            for gate_it in range(BNK_gates):
                truth_entry_bin=""
                for input_it in range(BNK_inputs):
                    gene_pos=gate_it*genes_per_gate+input_it

                    machine_state_pos=int(genes[gene_pos])
                    if machine_state_str[machine_state_pos] == "1":
                        truth_entry_bin+="1"
                    else:
                        truth_entry_bin+="0"

                    #Make sure we have the right amount of inputs for the truth table
                    if len(truth_entry_bin) > BNK_inputs:
                        raise Exception("Truth table entry "+truth_entry_bin.length()+" doesn't fit with the number of inputs "+BNK_inputs+".")

                truth_entry=int(truth_entry_bin, 2)
                truth_output_pos=gate_it*genes_per_gate+(BNK_inputs+truth_entry)
                next_state_str+=genes[truth_output_pos];
            next_state=int(next_state_str, 2)
            state_diagram.append(next_state)

        std_title="State transition diagram of individual #"+str(individual_line_num)
        std_file=destination_folder+"/"+std_title+".png"

        #Using Graphviz only
        graph=graphviz.Digraph(std_title, format='png')
        graph.attr(label=std_title)
        #graph.attr(fontsize='20')

        transition_it=0
        for transition in state_diagram:
            source=format(transition_it, format_rule)
            dest=format(transition, format_rule)
            graph.edge(source, dest)
            transition_it+=1

        graph.render(std_title)
        
        #should be working in a tmp dir, so move to official experiment folder
        copyfile(std_title+".png", std_file)
        os.remove(std_title)

    return


def PlotBNKLandscape(pop_folder,BNK_gates,BNK_inputs):
    pop_file = pop_folder+"/vre.csv"
    #destination_file=pop_folder+"/FitnessLandscape"+PLOT_EXTENSION
    destination_file_pgf=pop_folder+"/FitnessLandscape"+PGF_EXTENSION
    destination_file_png=pop_folder+"/FitnessLandscape"+PNG_EXTENSION

    X_AXIS="Symmetry"
    Y_AXIS="Robustness"


    #load data
    pop_data=pd.read_csv(pop_file)

    if set([X_AXIS,Y_AXIS]).issubset(pop_data.columns):
        fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

        ax.set_xlabel(X_AXIS)
        ax.set_ylabel(Y_AXIS)
        ax.set_zlabel('Fitness')

        ax.view_init(25,-135)

        X = pop_data[X_AXIS].values
        Y = pop_data[Y_AXIS].values
        Z = pop_data['Fitness'].values

        #Need to reduce dimensionality when BNK >= (4,4)
        Y2 = Y/np.max(Y)

        #ensure we have at least 3 distinct points on each axis
        if len(np.unique(X))>=3 and len(np.unique(Y))>=3:
            try:
                # Plot the surface.
                surf = ax.plot_trisurf(X, Y2, Z, cmap=cm.coolwarm, linewidth=0, antialiased=True)

            except Exception as e:
                print("Exception trying to plot the fitness landscape")
                print("Exception:"+str(e))
                print("File:"+pop_file)
                print("X:"+str(X))
                print("Y:"+str(Y))
                print("Z:"+str(Z))
        #else:
            #print("Not enough values to plot surface map")

        plt_title="BNK("+str(BNK_gates)+","+str(BNK_inputs)+") Fitness Landscape"
        plt.title(plt_title)

        plt.tight_layout()

        plt.savefig(destination_file_pgf)
        plt.savefig(destination_file_png)

        plt.close()
    return

def PlotFitnessLandscapeMultiModality(walk_folder,BNK_gates,BNK_inputs):

    flmm_title="BNK("+str(BNK_gates)+","+str(BNK_inputs)+") Fitness Landscape Multi-modality"
    #destination_file=walk_folder+"/"+flmm_title+PLOT_EXTENSION
    destination_file_pgf=walk_folder+"/"+flmm_title+PGF_EXTENSION
    destination_file_png=walk_folder+"/"+flmm_title+PNG_EXTENSION

    pop_file = walk_folder+"/vre.csv"
    #destination_file=pop_folder+"/FitnessLandscape"+PLOT_EXTENSION

    #load data
    pop_data=pd.read_csv(pop_file)

    #First, create the 2D matrix holding the grayscale image
    #initate every pixel to 0
    #iterate over all rows and update the corresponding pixel to the fitness value
    #plot the 2D grayscale image

    X_AXIS="Symmetry"
    Y_AXIS="Robustness"

    unique_x_values=pd.unique(pop_data[X_AXIS])
    unique_y_values=pd.unique(pop_data[Y_AXIS])

    unique_x_values.sort()
    unique_y_values.sort()

    pixel_x_bins=len(unique_x_values)
    pixel_y_bins=len(unique_y_values)

    #Max bins
    max_bins=10000
    x_steps=1
    y_steps=1

    if pixel_x_bins>max_bins:
        x_steps=math.floor(pixel_x_bins/max_bins)
        pixel_x_bins=max_bins

    if pixel_y_bins>max_bins:
        y_steps=math.floor(pixel_y_bins/max_bins)
        pixel_y_bins=max_bins

    sampled_unique_x_values=list()
    for i in range(pixel_x_bins):
        sampled_unique_x_values.append(unique_x_values[i*x_steps])


    sampled_unique_y_values=list()
    for i in range(pixel_y_bins):
        sampled_unique_y_values.append(unique_y_values[i*y_steps])

    #pixel_x_bins=max_bins if pixel_x_bins>max_bins else pixel_x_bins
    #pixel_y_bins=max_bins if pixel_y_bins>max_bins else pixel_y_bins

    peaks_img=np.zeros((pixel_x_bins, pixel_y_bins))

    for index, row in pop_data.iterrows():
        fitness=row["Fitness"]
        x_val=row[X_AXIS]
        y_val=row[Y_AXIS]

        try:
            #find X(Symmetry) bin
            #x_pos=math.floor(symmetry*pixel_x_bins)
            #x_pos = np.argmax(unique_x_values==x_val)
            x_pos=sampled_unique_x_values.index(x_val)

            #find Y(Robustness) bin
            #y_pos=math.floor(robustness*pixel_y_bins)
            #y_pos = np.argmax(unique_y_values==y_val)
            y_pos=sampled_unique_y_values.index(y_val)

            current_value=peaks_img[x_pos,y_pos]
            new_value=fitness if fitness>current_value else current_value
            peaks_img[x_pos,y_pos]=new_value

        except ValueError :
            pass


        #Make sure it fits in the subsampled matrix
        #if x_pos % x_steps==0 and y_pos % y_steps==0:
        #    sampled_x_pos=int(x_pos/x_steps)
        #    sampled_y_pos=int(y_pos/y_steps)

        #    current_value=peaks_img[sampled_x_pos,sampled_y_pos]
        #    new_value=fitness if fitness>current_value else current_value
        #    peaks_img[sampled_x_pos,sampled_y_pos]=new_value


    # the axis needs to be swapped (probably something to do with the custom matrix used to detect the peaks, investigate later...)
    #plt.xlabel(X_AXIS)
    #plt.ylabel(Y_AXIS)
    plt.xlabel(Y_AXIS)
    plt.ylabel(X_AXIS)

    detected_peaks = detect_peaks(peaks_img)
    number_of_peaks=detected_peaks.sum()

    #estimated number of peaks based on sampling
    original_x_size=len(unique_x_values)
    original_y_size=len(unique_y_values)

    total_pixels=original_x_size*original_y_size
    sampled_pixels=pixel_x_bins*pixel_y_bins

    sampled_ratio=sampled_pixels/total_pixels

    estimated_total_peaks=number_of_peaks/sampled_ratio

    #peaks_str=" ("+str(number_of_peaks)+" peaks)"
    peaks_str=" ("+str(int(estimated_total_peaks))+" estimated peaks)"

    plt.title(flmm_title+peaks_str)
    #dont comment this, it plots the image
    #plt.imshow(peaks_img, cmap=cm.PuRd)
    ax2=plt.gca()
    ax2.imshow(peaks_img, cmap=cm.PuRd, origin="lower")

#source:https://stackoverflow.com/questions/13583153/how-to-zoomed-a-portion-of-image-and-insert-in-the-same-plot-in-matplotlib

    ##Add zoomed subplot
    ## this is another inset axes over the main axes
    #a = plt.axes([0.2, 0.6, .2, .2], facecolor='y')
    ##plt.plot(t[:len(r)], r)

    #Zoom region
    zoom_size=50
    x_zoom_start=int(pixel_x_bins/1.75)
    x_zoom_stop=x_zoom_start+zoom_size
    y_zoom_start=int(pixel_y_bins/1.75)
    y_zoom_stop=y_zoom_start+zoom_size

    #zoomed_peaks_img=peaks_img[x_zoom_start:x_zoom_stop,y_zoom_start:y_zoom_stop]
    #zoomed_peaks_img=peaks_img[y_zoom_start:y_zoom_stop,x_zoom_start:x_zoom_stop]
    #plt.imshow(zoomed_peaks_img, cmap=cm.PuRd)
    #plt.title('Zoomed region')
    ##plt.xlim(0, 0.2)
    #plt.xticks([])
    #plt.yticks([])

    zoom_ratio=pixel_x_bins/(zoom_size*2)

    axins2 = zoomed_inset_axes(ax2, zoom=zoom_ratio, loc=3)
    #axins2.imshow(zoomed_peaks_img, cmap=cm.PuRd, origin="lower")
    axins2.imshow(peaks_img, cmap=cm.PuRd, origin="lower")
    #axins2.imshow(zoomed_peaks_img, cmap=cm.PuRd)

    # sub region of the original image
    #x1, x2, y1, y2 = 0, 50, 0, 50
    #axins2.set_xlim(x_zoom_start, x_zoom_stop)
    #axins2.set_ylim(y_zoom_start, y_zoom_stop)

    axins2.set_xlim(x_zoom_start, x_zoom_stop)
    axins2.set_ylim(y_zoom_start, y_zoom_stop)
    

    # fix the number of ticks on the inset axes
    #axins2.yaxis.get_major_locator().set_params(nbins=7)
    #axins2.xaxis.get_major_locator().set_params(nbins=7)
    axins2.tick_params(labelleft=False, labelbottom=False)
    axins2.set_xticks([])
    axins2.set_yticks([])

    # draw a bbox of the region of the inset axes in the parent axes and
    # connecting lines between the bbox and the inset axes area
    mark_inset(ax2, axins2, loc1=2, loc2=4, fc="none", ec="0.5")
    #use loc1 and loc2 to identify which corners to connect (1,2,3,4)


    plt.savefig(destination_file_pgf)
    plt.savefig(destination_file_png)
    plt.close()

    return


#source:https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array
def detect_peaks(image):
    """
    Takes an image and detect the peaks using the local maximum filter.
    Returns a boolean mask of the peaks (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """

    # define an 8-connected neighborhood
    neighborhood = generate_binary_structure(2,2)

    #apply the local maximum filter; all pixel of maximal value 
    #in their neighborhood are set to 1
    local_max = maximum_filter(image, footprint=neighborhood)==image
    #local_max is a mask that contains the peaks we are 
    #looking for, but also the background.
    #In order to isolate the peaks we must remove the background from the mask.

    #we create the mask of the background
    background = (image==0)

    #a little technicality: we must erode the background in order to 
    #successfully subtract it form local_max, otherwise a line will 
    #appear along the background border (artifact of the local maximum filter)
    eroded_background = binary_erosion(background, structure=neighborhood, border_value=1)

    #we obtain the final mask, containing only peaks, 
    #by removing the background from the local_max mask (xor operation)
    detected_peaks = local_max ^ eroded_background

    return detected_peaks

def PlotFitnessStats(walk_folder,BNK_gates,BNK_inputs):
    fd_title="BNK("+str(BNK_gates)+","+str(BNK_inputs)+") Fitness Distribution"
    #destination_file=walk_folder+"/"+fd_title+PLOT_EXTENSION
    destination_file_pgf=walk_folder+"/"+fd_title+PGF_EXTENSION
    destination_file_png=walk_folder+"/"+fd_title+PNG_EXTENSION

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

    plt.title(fd_title)

    plt.tight_layout()

    plt.savefig(destination_file_pgf)
    plt.savefig(destination_file_png)
    plt.close()

    return


def PlotVREWalk(df_walk_data, walk_folder):
    #destination_file=walk_folder+"/VRE_plot"+PLOT_EXTENSION
    destination_file_pgf=walk_folder+"/VRE_plot"+PGF_EXTENSION
    destination_file_png=walk_folder+"/VRE_plot"+PNG_EXTENSION

    xmin=df_walk_data["neighborhood_distance"].min()
    xmax=df_walk_data["neighborhood_distance"].max() +1
    
    fig = plt.figure()

    gs = fig.add_gridspec(3, 2)
    ax_rg = fig.add_subplot(gs[0, 0])
    ax_rg.set_xticks(np.arange(xmin, xmax))

    ax_eg = fig.add_subplot(gs[0, 1])
    ax_eg.set_xticks(np.arange(xmin, xmax))

    ax_rp = fig.add_subplot(gs[1, 0])
    ax_rp.set_xticks(np.arange(xmin, xmax))

    ax_size = fig.add_subplot(gs[1, 1])
    ax_size.set_xticks(np.arange(xmin, xmax))

    ax_nn = fig.add_subplot(gs[2, 0])
    ax_nn.set_xticks(np.arange(xmin, xmax))

    ax_up = fig.add_subplot(gs[2, 1])
    ax_up.set_xticks(np.arange(xmin, xmax))

    #Normalize some columns
    df_walk_data["total_neutral_neighbors_norm"] = df_walk_data.apply(lambda row: row["total_neutral_neighbors"] / row["total_neighborhood_size"], axis=1)
    df_walk_data["total_unique_phenotypes_norm"] = df_walk_data.apply(lambda row: row["total_unique_phenotypes"] / row["total_neighborhood_size"], axis=1)

    df_walk_data.plot(x="neighborhood_distance",y="local_genotype_robustness",ax=ax_rg,legend=None,title="Genotype robustness", xlabel="neighborhood", ylabel="genotype robustness")
    df_walk_data.plot(x="neighborhood_distance",y="local_genotype_evolvability",ax=ax_eg,legend=None,title="Genotype evolvability", xlabel="neighborhood", ylabel="genotype evolvability")
    df_walk_data.plot(x="neighborhood_distance",y="local_phenotype_robustness",ax=ax_rp,legend=None,title="Phenotype robustness", xlabel="neighborhood", ylabel="phenotype robustness")

    ax_nn_line=ax_nn.plot(df_walk_data["neighborhood_distance"],df_walk_data["total_neutral_neighbors"], label='absolute')
    ax_nn.set_title("Neutral neighbors")
    ax_nn.set_xlabel("neighborhood")
    ax_nn.set_ylabel("absolute")
    ax_nn_var = ax_nn.twinx()  # instantiate a second axes that shares the same x-axis
    ax_nn_var_line = ax_nn_var.plot(df_walk_data["neighborhood_distance"],df_walk_data["total_neutral_neighbors_norm"], label='relative', color='tab:red')
    ax_nn_var.set_ylabel("relative")

    # added these three lines
    ax_nn_lines = ax_nn_line + ax_nn_var_line
    labs = [line.get_label() for line in ax_nn_lines]
    ax_nn_var.legend(ax_nn_lines, labs)

    ax_up_line=ax_up.plot(df_walk_data["neighborhood_distance"],df_walk_data["total_unique_phenotypes"], label='absolute')
    ax_up.set_title("Unique phenotypes")
    ax_up.set_xlabel("neighborhood")
    ax_up.set_ylabel("absolute")
    ax_up_var = ax_up.twinx()  # instantiate a second axes that shares the same x-axis
    ax_up_var_line = ax_up_var.plot(df_walk_data["neighborhood_distance"],df_walk_data["total_unique_phenotypes_norm"], label='relative', color='tab:red')
    ax_up_var.set_ylabel("relative")

    # added these three lines
    ax_up_lines = ax_up_line + ax_up_var_line
    labs = [line.get_label() for line in ax_up_lines]
    ax_up_var.legend(ax_up_lines, labs)

    df_walk_data.plot(x="neighborhood_distance",y="total_neighborhood_size",ax=ax_size,legend=None,title="Neighborhood sizes", xlabel="neighborhood", ylabel="size")

    plt.tight_layout()

    plt.savefig(destination_file_pgf)
    plt.savefig(destination_file_png)
    plt.close()
    #https://timodenk.com/blog/exporting-matplotlib-plots-to-latex/

    return

def PlotVRESampleVariancePerNeighborhood(sample_folder):
    variance_files_prefix=sample_folder+"/vre_sample_variance_tracking_N*.csv"
    variance_files = glob.glob(variance_files_prefix)

    for file in variance_files:
        try:
            neighborhood = re.search('vre_sample_variance_tracking_N(.+?).csv', file).group(1)
        except AttributeError:
            neighborhood = '0' # apply your error handling

        #destination_file=sample_folder+"/VRE_plot_N"+str(neighborhood)+PLOT_EXTENSION
        destination_file_pgf=sample_folder+"/VRE_plot_N"+str(neighborhood)+PGF_EXTENSION
        destination_file_png=sample_folder+"/VRE_plot_N"+str(neighborhood)+PNG_EXTENSION

        #load data
        pop_variance_data=pd.read_csv(file)

        gs = gridspec.GridSpec(2, 2)
        ax_rg = plt.subplot(gs[0, 0])
        ax_eg = plt.subplot(gs[0, 1])
        ax_rp = plt.subplot(gs[1, 0])
        ax_ep = plt.subplot(gs[1, 1])
        
        var_color = 'tab:red'

        ax_rg.errorbar(pop_variance_data["count"], pop_variance_data["mean_genotype_robustness"], yerr=pop_variance_data["variance_genotype_robustness"], label='mean',  ecolor='silver')
        ax_rg.set_title("Genotype robustness")
        ax_rg.set_xlabel("Number of walks")
        ax_rg.set_ylabel("Mean")
        ax_rg.legend()

        ax_eg.errorbar(pop_variance_data["count"], pop_variance_data["mean_genotype_evolvability"], yerr=pop_variance_data["variance_genotype_evolvability"], label='mean',  ecolor='silver')
        ax_eg.set_title("Genotype evolvability")
        ax_eg.set_xlabel("Number of walks")
        ax_eg.set_ylabel("Mean")
        ax_eg.legend()

        ax_rp.errorbar(pop_variance_data["count"], pop_variance_data["mean_phenotype_robustness"], yerr=pop_variance_data["variance_phenotype_robustness"], label='mean',  ecolor='silver')
        ax_rp.set_title("Phenotype robustness")
        ax_rp.set_xlabel("Number of walks")
        ax_rp.set_ylabel("Mean")
        ax_rp.legend()

        ax_ep.errorbar(pop_variance_data["count"], pop_variance_data["mean_phenotype_evolvability"], yerr=pop_variance_data["variance_phenotype_evolvability"], label='mean',  ecolor='silver')
        ax_ep.set_title("Phenotype evolvability")
        ax_ep.set_xlabel("Number of walks")
        ax_ep.set_ylabel("Mean")
        ax_ep.legend()

        plt.tight_layout()
        plt.savefig(destination_file_pgf)
        plt.savefig(destination_file_png)
        plt.close()

    return

def PlotVRESampleVarianceAllNeighborhoods(sample_folder):
    variance_files_prefix=sample_folder+"/vre_sample_variance_tracking_N*.csv"
    variance_files = glob.glob(variance_files_prefix)

    variance_files.sort()

    #destination_file=sample_folder+"/VRE_plot_ALL"+PLOT_EXTENSION
    destination_file=sample_folder+"/VRE_plot_ALL"

    gs = gridspec.GridSpec(2, 2)
    
    ax_rg = plt.subplot(gs[0, 0])
    ax_rg.set_title("Genotype robustness")
    ax_rg.set_xlabel("Number of walks")
    ax_rg.set_ylabel("Mean")
    
    ax_eg = plt.subplot(gs[0, 1])
    ax_eg.set_title("Genotype evolvability")
    ax_eg.set_xlabel("Number of walks")
    ax_eg.set_ylabel("Mean")
    
    ax_rp = plt.subplot(gs[1, 0])
    ax_rp.set_title("Phenotype robustness")
    ax_rp.set_xlabel("Number of walks")
    ax_rp.set_ylabel("Mean")
    
    ax_ep = plt.subplot(gs[1, 1])
    ax_ep.set_title("Phenotype evolvability")
    ax_ep.set_xlabel("Number of walks")
    ax_ep.set_ylabel("Mean")
    

    for file in variance_files:
        try:
            neighborhood = re.search('vre_sample_variance_tracking_N(.+?).csv', file).group(1)
        except AttributeError:
            neighborhood = '0' # apply your error handling

        #load data
        pop_variance_data=pd.read_csv(file)

        plot_legend_label="N"+str(neighborhood)

        ax_rg.plot(pop_variance_data["count"], pop_variance_data["mean_genotype_robustness"], label=plot_legend_label)
        ax_eg.plot(pop_variance_data["count"], pop_variance_data["mean_genotype_evolvability"], label=plot_legend_label)
        ax_rp.plot(pop_variance_data["count"], pop_variance_data["mean_phenotype_robustness"], label=plot_legend_label)
        ax_ep.plot(pop_variance_data["count"], pop_variance_data["mean_phenotype_evolvability"], label=plot_legend_label)

    ax_rg_handles, ax_rg_labels = ax_rg.get_legend_handles_labels()
    ax_rg_labels, ax_rg_handles = zip(*sorted(zip(ax_rg_labels, ax_rg_handles), key=lambda t: t[0]))
    ax_rg.legend(ax_rg_handles, ax_rg_labels)

    ax_eg_handles, ax_eg_labels = ax_eg.get_legend_handles_labels()
    ax_eg_labels, ax_eg_handles = zip(*sorted(zip(ax_eg_labels, ax_eg_handles), key=lambda t: t[0]))
    ax_eg.legend(ax_eg_handles, ax_eg_labels)

    ax_rp_handles, ax_rp_labels = ax_rp.get_legend_handles_labels()
    ax_rp_labels, ax_rp_handles = zip(*sorted(zip(ax_rp_labels, ax_rp_handles), key=lambda t: t[0]))
    ax_rp.legend(ax_rp_handles, ax_rp_labels)

    ax_ep_handles, ax_ep_labels = ax_ep.get_legend_handles_labels()
    ax_ep_labels, ax_ep_handles = zip(*sorted(zip(ax_ep_labels, ax_ep_handles), key=lambda t: t[0]))
    ax_ep.legend(ax_ep_handles, ax_ep_labels)

    plt.tight_layout()
    plt.savefig(destination_file+PGF_EXTENSION)
    plt.savefig(destination_file+PNG_EXTENSION)
    plt.close()

    return

def PlotVREExperiment(df_experiment_data, experiment_folder):
    #destination_file=experiment_folder+"/VRE_plot_ALL"+PLOT_EXTENSION
    destination_file=experiment_folder+"/VRE_plot_ALL"

    gs = gridspec.GridSpec(2, 2)
    
    ax_rg = plt.subplot(gs[0, 0])
    ax_rg.set_title("Genotype robustness")
    ax_rg.set_xlabel("Number of samples")
    ax_rg.set_ylabel("Mean")
    ax_rg_y=defaultdict(list)#{Neighborhood:[values]}
    ax_rg_y_mean=defaultdict(list)
    
    ax_eg = plt.subplot(gs[0, 1])
    ax_eg.set_title("Genotype evolvability")
    ax_eg.set_xlabel("Number of samples")
    ax_eg.set_ylabel("Mean")
    ax_eg_y=defaultdict(list)
    ax_eg_y_mean=defaultdict(list)
    
    ax_rp = plt.subplot(gs[1, 0])
    ax_rp.set_title("Phenotype robustness")
    ax_rp.set_xlabel("Number of samples")
    ax_rp.set_ylabel("Mean")
    ax_rp_y=defaultdict(list)
    ax_rp_y_mean=defaultdict(list)
    
    ax_ep = plt.subplot(gs[1, 1])
    ax_ep.set_title("Phenotype evolvability")
    ax_ep.set_xlabel("Number of samples")
    ax_ep.set_ylabel("Mean")
    ax_ep_y=defaultdict(list)
    ax_ep_y_mean=defaultdict(list)

    for sample_id in df_experiment_data:
        for index, row in df_experiment_data[sample_id].iterrows():
            neighborhood=int(row['neighborhood_distance'])

            ax_rg_y[neighborhood].append(row['genotype_robustness_mean'])
            ax_rg_y_mean[neighborhood].append(np.mean(ax_rg_y[neighborhood]))

            ax_eg_y[neighborhood].append(row['genotype_evolvability_mean'])
            ax_eg_y_mean[neighborhood].append(np.mean(ax_eg_y[neighborhood]))

            ax_rp_y[neighborhood].append(row['phenotype_robustness_mean'])
            ax_rp_y_mean[neighborhood].append(np.mean(ax_rp_y[neighborhood]))
            
            ax_ep_y[neighborhood].append(row['phenotype_evolvability_mean'])
            ax_ep_y_mean[neighborhood].append(np.mean(ax_ep_y[neighborhood]))

    for neighborhood in ax_rg_y_mean:
        x_count=range(1,len(ax_rg_y[neighborhood])+1)
        plot_legend_label="N"+str(neighborhood)
        ax_rg.plot(x_count, ax_rg_y_mean[neighborhood], label=plot_legend_label)
        ax_eg.plot(x_count, ax_eg_y_mean[neighborhood], label=plot_legend_label)
        ax_rp.plot(x_count, ax_rp_y_mean[neighborhood], label=plot_legend_label)
        ax_ep.plot(x_count, ax_ep_y_mean[neighborhood], label=plot_legend_label)

    ax_rg_handles, ax_rg_labels = ax_rg.get_legend_handles_labels()
    ax_rg_labels, ax_rg_handles = zip(*sorted(zip(ax_rg_labels, ax_rg_handles), key=lambda t: t[0]))
    ax_rg.legend(ax_rg_handles, ax_rg_labels, loc = "upper right")

    ax_eg_handles, ax_eg_labels = ax_eg.get_legend_handles_labels()
    ax_eg_labels, ax_eg_handles = zip(*sorted(zip(ax_eg_labels, ax_eg_handles), key=lambda t: t[0]))
    ax_eg.legend(ax_eg_handles, ax_eg_labels, loc = "upper right")

    ax_rp_handles, ax_rp_labels = ax_rp.get_legend_handles_labels()
    ax_rp_labels, ax_rp_handles = zip(*sorted(zip(ax_rp_labels, ax_rp_handles), key=lambda t: t[0]))
    ax_rp.legend(ax_rp_handles, ax_rp_labels, loc = "upper right")

    ax_ep_handles, ax_ep_labels = ax_ep.get_legend_handles_labels()
    ax_ep_labels, ax_ep_handles = zip(*sorted(zip(ax_ep_labels, ax_ep_handles), key=lambda t: t[0]))
    ax_ep.legend(ax_ep_handles, ax_ep_labels, loc = "upper right")

    plt.tight_layout()
    plt.savefig(destination_file+PGF_EXTENSION)
    plt.savefig(destination_file+PNG_EXTENSION)
    plt.close()

    return

def PlotBRWalk(df_walk_data, walk_folder):
    #destination_file=walk_folder+"/BR_plot"+PLOT_EXTENSION
    destination_file=walk_folder+"/BR_plot"
    sample_file = walk_folder+"/br_distances_sample.csv"

    #load data
    pop_sample=pd.read_csv(sample_file)

    x_max=np.max(pop_sample["genotype_distance"])
    y_max=np.max(pop_sample["phenotype_distance"])

    try:
        PlotBRStep(pop_sample, df_walk_data, x_max, y_max, destination_file, False)
    except Exception as e:
        print("Failed to plot BR walk")
        print("Caught exception:" + str(e))
        print(traceback.format_exc())

    return


def PlotBRWalkVariance(walk_folder):
    #variance_destination_file=walk_folder+"/BR_variance_plot"+PLOT_EXTENSION
    variance_destination_file=walk_folder+"/BR_variance_plot"
    variance_file = walk_folder+"/br_variance_tracking.csv"

   #load data
    pop_variance_data=pd.read_csv(variance_file)

    pop_variance_data["genotype_distance_stddev"] = pop_variance_data.apply(lambda row: math.sqrt(row["genotype_distance_variance"]), axis=1)
    pop_variance_data["phenotype_distance_stddev"] = pop_variance_data.apply(lambda row: math.sqrt(row["phenotype_distance_variance"]), axis=1)
    pop_variance_data["eigen_value_0_stddev"] = pop_variance_data.apply(lambda row: math.sqrt(row["eigen_value_0_variance"]), axis=1)
    pop_variance_data["eigen_value_1_stddev"] = pop_variance_data.apply(lambda row: math.sqrt(row["eigen_value_1_variance"]), axis=1)
    pop_variance_data["ellipse_angle_stddev"] = pop_variance_data.apply(lambda row: math.sqrt(row["ellipse_angle_variance"]), axis=1)

    pop_variance_data["ellipse_angle_rad_mean"]=pop_variance_data.apply(lambda row: math.radians(row["ellipse_angle_mean"]), axis=1)
    #pop_variance_data["ellipse_angle_rad_stddev"]=pop_variance_data.apply(lambda row: math.sqrt(row["ellipse_angle_rad_mean"]), axis=1)

    gs = gridspec.GridSpec(3, 2)
    ax_genotype = plt.subplot(gs[0, 0])
    ax_phenotype = plt.subplot(gs[0, 1])
    ax_eigen_0 = plt.subplot(gs[1, 0])
    ax_eigen_1 = plt.subplot(gs[1, 1])

    ax_angle = plt.subplot(gs[2, 0])
    ax_angle_rad = plt.subplot(gs[2, 1])

    ax_genotype.errorbar(pop_variance_data["count"], pop_variance_data["genotype_distance_mean"], yerr=pop_variance_data["genotype_distance_stddev"], label='mean',  ecolor='silver')
    ax_genotype.set_title("Genotype distance")
    ax_genotype.set_xlabel("Number of steps")
    ax_genotype.set_ylabel("Mean")
    ax_genotype.legend()

    ax_phenotype.errorbar(pop_variance_data["count"], pop_variance_data["phenotype_distance_mean"], yerr=pop_variance_data["phenotype_distance_stddev"], label='mean',  ecolor='silver')
    ax_phenotype.set_title("Phenotype distance")
    ax_phenotype.set_xlabel("Number of steps")
    ax_phenotype.set_ylabel("Mean")
    ax_phenotype.legend()

    ax_eigen_0.errorbar(pop_variance_data["count"], pop_variance_data["eigen_value_0_mean"], yerr=pop_variance_data["eigen_value_0_stddev"], label='mean',  ecolor='silver')
    ax_eigen_0.set_title("Eigen value 0")
    ax_eigen_0.set_xlabel("Number of steps")
    ax_eigen_0.set_ylabel("Mean")
    ax_eigen_0.legend()

    ax_eigen_1.errorbar(pop_variance_data["count"], pop_variance_data["eigen_value_1_mean"], yerr=pop_variance_data["eigen_value_1_stddev"], label='mean',  ecolor='silver')    
    ax_eigen_1.set_title("Eigen value 1")
    ax_eigen_1.set_xlabel("Number of steps")
    ax_eigen_1.set_ylabel("Mean")
    ax_eigen_1.legend()

    ax_angle.errorbar(pop_variance_data["count"], pop_variance_data["ellipse_angle_mean"], yerr=pop_variance_data["ellipse_angle_stddev"], label='mean',  ecolor='silver')    
    ax_angle.set_title("Ellipse angle in degree")
    ax_angle.set_xlabel("Number of steps")
    ax_angle.set_ylabel("Mean")
    ax_angle.legend()

    ax_angle_rad.plot(pop_variance_data["count"], pop_variance_data["ellipse_angle_rad_mean"], label='mean')    
    ax_angle_rad.set_title("Ellipse angle in radian")
    ax_angle_rad.set_xlabel("Number of steps")
    ax_angle_rad.set_ylabel("Mean")
    ax_angle_rad.legend()

    fig = plt.gcf()
    fig.set_size_inches(10, 8)
    plt.tight_layout()  # otherwise the right y-label is slightly clipped

    plt.savefig(variance_destination_file+PGF_EXTENSION)
    plt.savefig(variance_destination_file+PNG_EXTENSION)
    plt.close()
   
    return

def PlotBRWalkEllipseEvolution(walk_folder, steps=1):
    # define the name of the directory to be created
    ellipse_folder = walk_folder+"/ellipse"

    try:
        os.mkdir(ellipse_folder)
    except OSError:
        print ("Creation of the directory %s failed" % ellipse_folder)

    sample_file = walk_folder+"/br_distances_sample.csv"
    ellipse_file = walk_folder+"/br_ellipse_tracking.csv"

    pop_ellipse_data=pd.read_csv(ellipse_file)
    pop_sample=pd.read_csv(sample_file)

    x_max=np.max(pop_sample["genotype_distance"])
    y_max=np.max(pop_sample["phenotype_distance"])

    #Determine the intervals at which we need to generate an image
    #ellipse_step_size = math.floor(len(pop_ellipse_data.index) / steps)
    sample_step_size = math.floor(len(pop_sample.index) / steps)

    gif_file=ellipse_folder+"/BR_plot_movie.gif"
    gif_writer = imageio.get_writer(gif_file, mode='I')

    for step in np.arange(steps):
        sub_pop_ellipse_data=pop_ellipse_data.iloc[int((step+1)*sample_step_size)]
        sub_pop_sample=pop_sample.iloc[0:int((step+1)*sample_step_size)]
        
        #destination_file=ellipse_folder+"/BR_plot_step_"+str(step)+PLOT_EXTENSION
        destination_file=ellipse_folder+"/BR_plot_step_"+str(step)

        PlotBRStep(sub_pop_sample, sub_pop_ellipse_data, x_max, y_max, destination_file)

        image = imageio.imread(destination_file+PNG_EXTENSION)
        gif_writer.append_data(image)

        #Add the last image a few times to give some delay before restarting the loop
        if (step+1 >= steps):
            for x in range(0, 10):
                gif_writer.append_data(image)

        #Clean up
        os.remove(destination_file+PGF_EXTENSION)
        os.remove(destination_file+PNG_EXTENSION)

    plt.close()

    return

def PlotBRStep(sub_pop_sample, sub_pop_ellipse_data, x_max, y_max, destination_file, show_step=True):
    ax_br=plt.gca()

    # num of steps included
    total_steps=len(sub_pop_sample)

    # Create 2D Histogram plot
    #ax_br.hist2d(sub_pop_sample["genotype_distance"], sub_pop_sample["phenotype_distance"], bins=[np.arange(x_max+1) - 0.5, np.arange(y_max+1) - 0.5], alpha=0.9, cmap='Greys')
    if y_max>1.0:
        ax_br.hist2d(sub_pop_sample["genotype_distance"], sub_pop_sample["phenotype_distance"], bins=[np.arange(x_max+1) - 0.5, np.arange(y_max+1) - 0.5], alpha=0.9, cmap='Greys')
    else:
        ax_br.hist2d(sub_pop_sample["genotype_distance"], sub_pop_sample["phenotype_distance"], bins=[np.arange(x_max+1) - 0.5, np.arange(0,y_max+0.2,0.1)-0.05], alpha=0.9, cmap='Greys')

    sub_pop_sample.plot.scatter(x="genotype_distance",y="phenotype_distance",ax=ax_br, s=1, c='green')

    try:
        # Add confidence ellipses
                                                                        #Need to escape the percent sign if using LaTeX
        confidence_ellipse_from_stats(sub_pop_ellipse_data, ax_br, confidence=0.65, label='65\%', edgecolor=(1, 0, 0.2))
        confidence_ellipse_from_stats(sub_pop_ellipse_data, ax_br, confidence=0.90, label='90\%', edgecolor=(0.8, 0, 0.5))
        confidence_ellipse_from_stats(sub_pop_ellipse_data, ax_br, confidence=0.95, label='95\%', edgecolor=(0.5, 0, 0.9))
    except Exception as e:
        print("Failed to plot the confidence ellipse")
        print("Caught exception:" + str(e))
        print(traceback.format_exc())

    # Add eigen vectors
    eigen_vectors_from_stats(sub_pop_ellipse_data, ax_br)

    # Plot mean value
    ax_br.scatter(sub_pop_ellipse_data["combinations_genotype_distance_mean"], sub_pop_ellipse_data["combinations_phenotype_distance_mean"], c='red', s=5)

    if show_step:
        # Show current step
        current_step="Step: {}".format(total_steps)
        bbox_props = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.9)
        plt.text(0.98, 0.02, current_step, horizontalalignment='right', verticalalignment='bottom', transform=ax_br.transAxes, bbox=bbox_props)

    # Layout and titles
    ax_br.set_title("Baseline Robustness of a single walk")
    ax_br.set_xlabel("genotype distance")
    ax_br.set_ylabel("phenotype distance")
    ax_br.set_xlim([-1, x_max+1])
    ax_br.set_ylim([-1, y_max+1])
    #ax_br.legend(loc="upper left")

    plt.legend(loc="upper left")

    plt.savefig(destination_file+PGF_EXTENSION)
    plt.savefig(destination_file+PNG_EXTENSION)
    plt.close()

    return

def PlotBRSample(df_sample_data, sample_folder, track_variance):
    #destination_file=sample_folder+"/BR_plot"+PLOT_EXTENSION
    destination_file=sample_folder+"/BR_plot"

    gs = gridspec.GridSpec(3, 2)

    ax_genotype = plt.subplot(gs[0, 0])
    ax_genotype_y=[]
    ax_genotype_y_mean=[]
    ax_genotype_y_variance=[]

    ax_phenotype = plt.subplot(gs[0, 1])
    ax_phenotype_y=[]
    ax_phenotype_y_mean=[]
    ax_phenotype_y_variance=[]

    ax_eigen_0 = plt.subplot(gs[1, 0])
    ax_eigen_0_y=[]
    ax_eigen_0_y_mean=[]
    ax_eigen_0_y_variance=[]

    ax_eigen_1 = plt.subplot(gs[1, 1])
    ax_eigen_1_y=[]
    ax_eigen_1_y_mean=[]
    ax_eigen_1_y_variance=[]

    ax_angle = plt.subplot(gs[2, 0])
    ax_angle.set_ylim([-180, 180])
    ax_angle_y=[]
    ax_angle_y_mean=[]
    ax_angle_y_variance=[]

    ax_angle_rad = plt.subplot(gs[2, 1])
    ax_angle_rad_y=[]
    ax_angle_rad_y_mean=[]
    ax_angle_rad_y_variance=[]

    #for pop in populations:
    for walk_id in df_sample_data:
        ax_genotype_y.append(df_sample_data[walk_id]["combinations_genotype_distance_mean"])
        ax_genotype_y_mean.append(np.mean(ax_genotype_y))
        ax_genotype_y_variance.append(np.var(ax_genotype_y))

        ax_phenotype_y.append(df_sample_data[walk_id]["combinations_phenotype_distance_mean"])
        ax_phenotype_y_mean.append(np.mean(ax_phenotype_y))
        ax_phenotype_y_variance.append(np.var(ax_phenotype_y))

        ax_eigen_0_y.append(df_sample_data[walk_id]["combinations_eigen_value_0"])
        ax_eigen_0_y_mean.append(np.mean(ax_eigen_0_y))
        ax_eigen_0_y_variance.append(np.var(ax_eigen_0_y))

        ax_eigen_1_y.append(df_sample_data[walk_id]["combinations_eigen_value_1"])
        ax_eigen_1_y_mean.append(np.mean(ax_eigen_1_y))
        ax_eigen_1_y_variance.append(np.var(ax_eigen_1_y))

        ax_angle_y.append(df_sample_data[walk_id]["ellipse_angle"])
        ax_angle_y_mean.append(np.mean(ax_angle_y))
        ax_angle_y_variance.append(np.var(ax_angle_y))

        ax_angle_rad_y.append(math.radians(df_sample_data[walk_id]["ellipse_angle"]))
        ax_angle_rad_y_mean.append(np.mean(ax_angle_rad_y))
        ax_angle_rad_y_variance.append(np.var(ax_angle_rad_y))

    x_count=range(1,len(ax_genotype_y)+1)
    
    ax_genotype.errorbar(x_count, ax_genotype_y_mean, yerr=ax_genotype_y_variance, label='mean',  ecolor='silver')
    ax_genotype.set_title("Genotype distance")
    ax_genotype.set_xlabel("Number of walks")
    ax_genotype.set_ylabel("Mean")
    ax_genotype.legend()

    ax_phenotype.errorbar(x_count, ax_phenotype_y_mean, yerr=ax_phenotype_y_variance, label='mean',  ecolor='silver')
    ax_phenotype.set_title("Phenotype distance")
    ax_phenotype.set_xlabel("Number of walks")
    ax_phenotype.set_ylabel("Mean")
    ax_phenotype.legend()

    ax_eigen_0.errorbar(x_count, ax_eigen_0_y_mean, yerr=ax_eigen_0_y_variance, label='mean',  ecolor='silver')
    ax_eigen_0.set_title("Eigen value 0")
    ax_eigen_0.set_xlabel("Number of walks")
    ax_eigen_0.set_ylabel("Mean")
    ax_eigen_0.legend()

    ax_eigen_1.errorbar(x_count, ax_eigen_1_y_mean, yerr=ax_eigen_1_y_variance, label='mean',  ecolor='silver')    
    ax_eigen_1.set_title("Eigen value 1")
    ax_eigen_1.set_xlabel("Number of walks")
    ax_eigen_1.set_ylabel("Mean")
    ax_eigen_1.legend()

    ax_angle.errorbar(x_count, ax_angle_y_mean, yerr=ax_angle_y_variance, label='mean',  ecolor='silver')    
    ax_angle.set_title("Ellipse angle in degree")
    ax_angle.set_xlabel("Number of walks")
    ax_angle.set_ylabel("Mean")
    ax_angle.legend()

    ax_angle_rad.errorbar(x_count, ax_angle_rad_y_mean, yerr=ax_angle_rad_y_variance, label='mean',  ecolor='silver')    
    ax_angle_rad.set_title("Ellipse angle in radian")
    ax_angle_rad.set_xlabel("Number of walks")
    ax_angle_rad.set_ylabel("Mean")
    ax_angle_rad.legend()

    fig = plt.gcf()
    fig.set_size_inches(10, 8)
    plt.tight_layout()  # otherwise the right y-label is slightly clipped

    plt.savefig(destination_file+PGF_EXTENSION)
    plt.savefig(destination_file+PNG_EXTENSION)
    plt.close()
    #https://timodenk.com/blog/exporting-matplotlib-plots-to-latex/
    return

def PlotBRExperiment(df_experiment_data, experiment_folder):
    #destination_file=experiment_folder+"/BR_plot"+PLOT_EXTENSION
    destination_file=experiment_folder+"/BR_plot"

    gs = gridspec.GridSpec(3, 2)

    ax_genotype = plt.subplot(gs[0, 0])
    ax_genotype_y=[]
    ax_genotype_y_mean=[]
    ax_genotype_y_variance=[]

    ax_phenotype = plt.subplot(gs[0, 1])
    ax_phenotype_y=[]
    ax_phenotype_y_mean=[]
    ax_phenotype_y_variance=[]

    ax_eigen_0 = plt.subplot(gs[1, 0])
    ax_eigen_0_y=[]
    ax_eigen_0_y_mean=[]
    ax_eigen_0_y_variance=[]

    ax_eigen_1 = plt.subplot(gs[1, 1])
    ax_eigen_1_y=[]
    ax_eigen_1_y_mean=[]
    ax_eigen_1_y_variance=[]

    ax_angle = plt.subplot(gs[2, 0])
    ax_angle.set_ylim([-180, 180])
    ax_angle_y=[]
    ax_angle_y_mean=[]
    ax_angle_y_variance=[]

    ax_angle_rad = plt.subplot(gs[2, 1])
    ax_angle_rad_y=[]
    ax_angle_rad_y_mean=[]
    ax_angle_rad_y_variance=[]

    for sample_id in df_experiment_data:
        ax_genotype_y.append(df_experiment_data[sample_id]["genotype_distance_mean"])
        ax_genotype_y_mean.append(np.mean(ax_genotype_y))
        ax_genotype_y_variance.append(np.var(ax_genotype_y))

        ax_phenotype_y.append(df_experiment_data[sample_id]["phenotype_distance_mean"])
        ax_phenotype_y_mean.append(np.mean(ax_phenotype_y))
        ax_phenotype_y_variance.append(np.var(ax_phenotype_y))

        ax_eigen_0_y.append(df_experiment_data[sample_id]["eigen_value_0_mean"])
        ax_eigen_0_y_mean.append(np.mean(ax_eigen_0_y))
        ax_eigen_0_y_variance.append(np.var(ax_eigen_0_y))

        ax_eigen_1_y.append(df_experiment_data[sample_id]["eigen_value_1_mean"])
        ax_eigen_1_y_mean.append(np.mean(ax_eigen_1_y))
        ax_eigen_1_y_variance.append(np.var(ax_eigen_1_y))

        ax_angle_y.append(df_experiment_data[sample_id]["ellipse_angle_mean"])
        ax_angle_y_mean.append(np.mean(ax_angle_y))
        ax_angle_y_variance.append(np.var(ax_angle_y))

        ax_angle_rad_y.append(math.radians(df_experiment_data[sample_id]["ellipse_angle_mean"]))
        ax_angle_rad_y_mean.append(np.mean(ax_angle_rad_y))
        ax_angle_rad_y_variance.append(np.var(ax_angle_rad_y))


    x_count=range(1,len(ax_genotype_y)+1)

    ax_genotype.errorbar(x_count, ax_genotype_y_mean, yerr=ax_genotype_y_variance, label='mean',  ecolor='silver')
    ax_genotype.set_title("Genotype distance")
    ax_genotype.set_xlabel("Number of samples")
    ax_genotype.set_ylabel("Mean")
    ax_genotype.legend()

    ax_phenotype.errorbar(x_count, ax_phenotype_y_mean, yerr=ax_phenotype_y_variance, label='mean',  ecolor='silver')
    ax_phenotype.set_title("Phenotype distance")
    ax_phenotype.set_xlabel("Number of samples")
    ax_phenotype.set_ylabel("Mean")
    ax_phenotype.legend()

    ax_eigen_0.errorbar(x_count, ax_eigen_0_y_mean, yerr=ax_eigen_0_y_variance, label='mean',  ecolor='silver')
    ax_eigen_0.set_title("Eigen value 0")
    ax_eigen_0.set_xlabel("Number of samples")
    ax_eigen_0.set_ylabel("Mean")
    ax_eigen_0.legend()

    ax_eigen_1.errorbar(x_count, ax_eigen_1_y_mean, yerr=ax_eigen_1_y_variance, label='mean',  ecolor='silver')    
    ax_eigen_1.set_title("Eigen value 1")
    ax_eigen_1.set_xlabel("Number of samples")
    ax_eigen_1.set_ylabel("Mean")
    ax_eigen_1.legend()

    ax_angle.errorbar(x_count, ax_angle_y_mean, yerr=ax_angle_y_variance, label='mean',  ecolor='silver')    
    ax_angle.set_title("Ellipse angle in degree")
    ax_angle.set_xlabel("Number of samples")
    ax_angle.set_ylabel("Mean")
    ax_angle.legend()

    ax_angle_rad.errorbar(x_count, ax_angle_rad_y_mean, yerr=ax_angle_rad_y_variance, label='mean',  ecolor='silver')    
    ax_angle_rad.set_title("Ellipse angle in radian")
    ax_angle_rad.set_xlabel("Number of samples")
    ax_angle_rad.set_ylabel("Mean")
    ax_angle_rad.legend()

    fig = plt.gcf()
    fig.set_size_inches(10, 8)
    plt.tight_layout()  # otherwise the right y-label is slightly clipped

    plt.savefig(destination_file+PGF_EXTENSION)
    plt.savefig(destination_file+PNG_EXTENSION)
    plt.close()
    #https://timodenk.com/blog/exporting-matplotlib-plots-to-latex/
    return

def PlotBEWalkVariance(walk_folder):
    #variance_destination_file=walk_folder+"/BE_variance_plot"+PLOT_EXTENSION
    variance_destination_file=walk_folder+"/BE_variance_plot"
    variance_file = walk_folder+"/be_variance_tracking.csv"

    pop_variance_data=pd.read_csv(variance_file)

    ax_be = plt.gca()

    err_handles = ax_be.errorbar(pop_variance_data["count"], pop_variance_data["mean"], yerr=pop_variance_data["variance"], label='mean',  ecolor='silver')
    ax_be.set_title("Baseline evolvability for a single walk")
    ax_be.set_xlabel("Number of steps")
    ax_be.set_ylabel("Mean")

    ax_be_closeness = ax_be.twinx()  # instantiate a second axes that shares the same x-axis
    ax_be_closeness_color = 'tab:green'
    line3 = ax_be_closeness.plot(pop_variance_data["count"], pop_variance_data["closeness"], label='closeness', color=ax_be_closeness_color)
    ax_be_closeness.set_ylabel("Closeness")

    # added these three lines
    # get handles
    handles, labels = ax_be.get_legend_handles_labels()
    # remove the errorbars
    lines = [h[0] for h in handles]
    lines = handles+line3
    labs = [line.get_label() for line in lines]
    # use them in the legend
    plt.legend(lines, labs, loc='upper right')

    plt.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.savefig(variance_destination_file+PGF_EXTENSION)
    plt.savefig(variance_destination_file+PNG_EXTENSION)
    plt.close()
   
    return

#source:https://stackoverflow.com/questions/8500700/how-to-plot-a-gradient-color-line-in-matplotlib\
def colorline(ax, x, y, z=None, cmap=plt.get_cmap('copper'), norm=plt.Normalize(0.0, 1.0), linewidth=3, alpha=1.0):
    """
    http://nbviewer.ipython.org/github/dpsanders/matplotlib-examples/blob/master/colorline.ipynb
    http://matplotlib.org/examples/pylab_examples/multicolored_line.html
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    """

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))

    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])

    z = np.asarray(z)

    segments = make_segments(x, y)
    lc = mcoll.LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)

    #ax = plt.gca()
    ax.add_collection(lc)

    return lc


def make_segments(x, y):
    """
    Create list of line segments from x and y coordinates, in the correct format
    for LineCollection: an array of the form numlines x (points per line) x 2 (x
    and y) array
    """

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments


#BNK(3,2) example:
# -S 0.375 -R 0.8660254037844387
def PlotBNK_BE(walk_folder,target_sym=0.375,target_rob=0.8660254037844387):

    pop_file = walk_folder+"/vre.csv"

    #load data
    pop_data=pd.read_csv(pop_file)

    #Skip if too many points... (quick hack to discard big experiements)
    if len(pop_data)>500:
        return

    #Get list of genotypes traversed
    unique_walk_data = pop_data.drop_duplicates(subset=['Genotype'])

    title_bnk="Baseline evolvability"
    #destination_file=walk_folder+"/BNK_BE"+PLOT_EXTENSION
    destination_file=walk_folder+"/BNK_BE"
    ax_bnk=plt.gca()

    plt.xlim(-0.05,1.05)
    plt.ylim(0.5,1.05)

    plt.xlabel("Symmetry")
    plt.ylabel("Robustness")
    plt.title(title_bnk)


    default_color="black"
    start_color="blue"
    end_color="red"
    target_color="gold"

    #Setup target point
    #ax_bnk.plot(target_sym,target_rob,marker="*",markersize=12,color=end_color)

    #Mark starting point
    ax_bnk.plot(unique_walk_data["Symmetry"].iloc[0],unique_walk_data["Robustness"].iloc[0],marker="s",markersize=6,color=start_color)
    #ax_bnk.plot(unique_walk_data["Symmetry"][0],unique_walk_data["Robustness"][0],marker="D",markersize=15)

    #Mark stopping point
    last_point=len(unique_walk_data["Symmetry"])-1
    last_sym=unique_walk_data["Symmetry"].iloc[last_point]
    last_rob=unique_walk_data["Robustness"].iloc[last_point]
    if(last_sym != target_sym or last_rob != target_rob):
        #ax_bnk.plot(last_sym,last_rob,marker=6,markersize=15)
        ax_bnk.plot(last_sym,last_rob,marker="D",markersize=6,color=end_color)


    ax_bnk.scatter(unique_walk_data["Symmetry"],unique_walk_data["Robustness"],color=default_color,s=6)

    #Custom line plot with changing line colors
    x=unique_walk_data["Symmetry"]
    y=unique_walk_data["Robustness"]
    path = mpath.Path(np.column_stack([x, y]))
    verts = path.interpolated(steps=3).vertices
    x, y = verts[:, 0], verts[:, 1]
    z = np.linspace(0, 1, len(x))
    lines=colorline(ax_bnk, x, y, z, cmap=plt.get_cmap('coolwarm'), linewidth=1)

    ######### ARROWS
    #line=ax_bnk.plot(unique_walk_data["Symmetry"],unique_walk_data["Robustness"],color=default_color)[0]
    #Add arrow direction 
    #if(len(unique_walk_data)>1):
    #    arrow_color = lines.get_color()
    #    for d in range(len(unique_walk_data)-1):
    #        line.axes.annotate('', xytext=(unique_walk_data["Symmetry"].iloc[d], unique_walk_data["Robustness"].iloc[d]), xy=(unique_walk_data["Symmetry"].iloc[d+1], unique_walk_data["Robustness"].iloc[d+1]), arrowprops=dict(arrowstyle="->", color=arrow_color), size=15)


    ########## EDGE INDEX
    ##Prepare point annotation
    #points_annotation=defaultdict(list)
    #unique_walk_data["key"] = unique_walk_data['Symmetry'].astype(str) +";"+ unique_walk_data["Robustness"].astype(str)
    #for d in range(len(unique_walk_data)-1):
    #    point=unique_walk_data["key"].iloc[d]
    #    points_annotation[point].append(str(d))
    #for key, value in points_annotation.items():
    #    #txt_value=','.join(value)
    #    txt_value=value[0]
    #    x,y=key.split(";")
    #    plt.annotate(txt_value, (float(x), float(y)))


    #Prepare point annotation
    points_annotation=defaultdict(list)
    #unique_walk_data["key"] = unique_walk_data['Symmetry'].astype(str) +";"+ unique_walk_data["Robustness"].astype(str)
    #Change to fix SettingWithCopyWarning
    fn = lambda row: str(row["Symmetry"]) +";"+ str(row["Robustness"]) # define a function for the new column
    col = unique_walk_data.apply(fn, axis=1) # get column data with an index
    unique_walk_data = unique_walk_data.assign(key=col.values) # assign values to column 'c'

    #points=unique_walk_data["key"].unique().tolist()

    for d in range(len(unique_walk_data)):
        #print(d)
        point=unique_walk_data["key"].iloc[d]
        points_annotation[point].append(str(d))

    visited_vertex=set()
    step=0

    bgcolor=(1, 1, 1, 0.35)

    for key, value in points_annotation.items():
        if key not in visited_vertex:
            #txt_value=','.join(value)
            #txt_value=value[0]
            x,y=key.split(";")        
            visited_vertex.add(key)
            #plt.annotate(txt_value, (float(x), float(y)))
            plt.annotate(str(step), (float(x)+0.005, float(y)+0.005),fontsize=8, backgroundcolor=bgcolor)
            step=step+1

    
    plt.tight_layout()

    plt.savefig(destination_file+PGF_EXTENSION)
    plt.savefig(destination_file+PNG_EXTENSION)
    plt.close()

    return

def PlotBESample(df_sample_data, sample_folder, track_variance):
    #destination_file=sample_folder+"/BE_plot"+PLOT_EXTENSION
    destination_file=sample_folder+"/BE_plot"

    y=[]
    y_mean=[]
    y_variance=[]

    for walk_id in df_sample_data:
        y.append(df_sample_data[walk_id]["baseline_value"])
        y_mean.append(np.mean(y))
        y_variance.append(np.var(y))

    ax_be = plt.gca()
    
    x_count=range(1,len(y)+1)

    ax_be.plot(x_count,y, label='value')
    if(track_variance):
        ax_be.errorbar(x_count,y_mean, yerr=y_variance, label='mean',  ecolor='silver')
    ax_be.legend()

    # Layout and titles
    ax_be.set_title("Evolution of baseline evolvability over multiple walks")
    ax_be.set_xlabel("Number of walks")
    ax_be.set_ylabel("Baseline evolvability")

    plt.savefig(destination_file+PGF_EXTENSION)
    plt.savefig(destination_file+PNG_EXTENSION)
    plt.close()
    #https://timodenk.com/blog/exporting-matplotlib-plots-to-latex/
    return

def PlotBEExperiment(df_experiment_data, experiment_folder, track_variance):
    #destination_file=experiment_folder+"/BE_plot"+PLOT_EXTENSION
    destination_file=experiment_folder+"/BE_plot"

    y=[]
    y_mean=[]
    y_variance=[]

    for sample_id in df_experiment_data:
        y.append(df_experiment_data[sample_id]["baseline_evolvability_mean"])
        y_mean.append(np.mean(y))
        y_variance.append(np.var(y))

    ax_be = plt.gca()
    
    x_count=range(1,len(y)+1)

    ax_be.plot(x_count,y, label='value')
    if(track_variance):
        ax_be.errorbar(x_count,y_mean, yerr=y_variance, label='mean',  ecolor='silver')
    ax_be.legend()

    # Layout and titles
    ax_be.set_title("Evolution of baseline evolvability over multiple samples")
    ax_be.set_xlabel("Number of samples")
    ax_be.set_ylabel("Baseline evolvability")

    plt.savefig(destination_file+PGF_EXTENSION)
    plt.savefig(destination_file+PNG_EXTENSION)
    plt.close()
    #https://timodenk.com/blog/exporting-matplotlib-plots-to-latex/
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


##source:https://matplotlib.org/devdocs/gallery/statistics/confidence_ellipse.html
#def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
#    if x.size != y.size:
#        raise ValueError("x and y must be the same size")

#    cov_mat = np.cov(x, y)
#    pearson = cov_mat[0, 1]/np.sqrt(cov_mat[0, 0] * cov_mat[1, 1])
#    # Using a special case to obtain the eigenvalues of this two-dimensionl dataset.
#    ell_radius_x = np.sqrt(1 + pearson)
#    ell_radius_y = np.sqrt(1 - pearson)
#    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, facecolor=facecolor, **kwargs)

#    # Calculating the stdandard deviation of x from the squareroot of the variance and multiplying with the given number of standard deviations.
#    scale_x = np.sqrt(cov_mat[0, 0]) * n_std
#    mean_x = np.mean(x)

#    # calculating the stdandard deviation of y ...
#    scale_y = np.sqrt(cov_mat[1, 1]) * n_std
#    mean_y = np.mean(y)

#    transf = transforms.Affine2D() \
#        .rotate_deg(45) \
#        .scale(scale_x, scale_y) \
#        .translate(mean_x, mean_y)

#    ellipse.set_transform(transf + ax.transData)
#    return ax.add_patch(ellipse)

def confidence_ellipse_from_stats(pop_stats, ax, confidence=0.95, facecolor='none', **kwargs):
    cov00=float(pop_stats["combinations_covariance_0_0"])
    cov01=float(pop_stats["combinations_covariance_0_1"])
    cov10=float(pop_stats["combinations_covariance_1_0"])
    cov11=float(pop_stats["combinations_covariance_1_1"])
    cov_mat = np.array([[cov00, cov01], [cov10, cov11]], np.float64)
    pearson = cov_mat[0, 1]/np.sqrt(cov_mat[0, 0] * cov_mat[1, 1])

    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    
    if math.isnan(ell_radius_y) or math.isnan(ell_radius_x):
        raise ValueError("Input is not a positive integer")
    

    ellipse = Ellipse((0, 0), width=ell_radius_x * 2, height=ell_radius_y * 2, facecolor=facecolor, **kwargs)

    chisquare_val= chi2.ppf(confidence,2)

    scale_x = np.sqrt(cov_mat[0, 0] * chisquare_val)
    mean_x = float(pop_stats["combinations_genotype_distance_mean"])

    scale_y = np.sqrt(cov_mat[1, 1] * chisquare_val)
    mean_y = float(pop_stats["combinations_phenotype_distance_mean"])

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)

    ax.add_patch(ellipse)

    return

def eigen_vectors_from_stats(pop_stats, ax, **kwargs):
    eigen_vec_00=float(pop_stats["combinations_eigen_vector_0_0"])
    eigen_vec_01=float(pop_stats["combinations_eigen_vector_0_1"])
    eigen_vec_10=float(pop_stats["combinations_eigen_vector_1_0"])
    eigen_vec_11=float(pop_stats["combinations_eigen_vector_1_1"])
    eig_vecs = np.array([[eigen_vec_00, eigen_vec_01],[eigen_vec_10, eigen_vec_11]])

    eigen_val_0 = float(pop_stats["combinations_eigen_value_0"])
    eigen_val_1 = float(pop_stats["combinations_eigen_value_1"])
    eig_vals = np.array([eigen_val_0, eigen_val_1])

    mean_x = float(pop_stats["combinations_genotype_distance_mean"])
    mean_y = float(pop_stats["combinations_phenotype_distance_mean"])

    # Get the largest eigenvalue
    largest_eigenval = max(eig_vals)

    # Get the index of the largest eigenvector
    largest_eigenvec_ind_c = np.argwhere(eig_vals == max(eig_vals))[0][0]
    largest_eigenvec = eig_vecs[:,largest_eigenvec_ind_c]


    # Get the smallest eigenvector and eigenvalue
    smallest_eigenval = min(eig_vals)
    if largest_eigenvec_ind_c == 0:
        smallest_eigenvec = eig_vecs[:,1]
    else:
        smallest_eigenvec = eig_vecs[:,0]


    #soa = np.array([[mean_x, mean_y, 
    #                 eig_vals[0] * eig_vecs[0][0], 
    #                 eig_vals[0] * eig_vecs[1][0]]])

    #soa1 = np.array([[mean_x, mean_y, 
    #                  eig_vals[1] * eig_vecs[0][1], 
    #                  eig_vals[1] * eig_vecs[1][1]]])

    #X, Y, U, V = zip(*soa)
    #X1, Y1, U1, V1 = zip(*soa1)

    #quiveropts = dict(headaxislength=0, color='red', headlength=0, units='xy', scale=6)
    quiveropts = dict(headaxislength=0, color='red', headlength=0, units='xy',angles='xy',scale=1)

    #plt.quiver(X, Y, U, V, **quiveropts)
    #plt.quiver(X1, Y1, U1, V1, **quiveropts)

    plt.quiver(mean_x, mean_y, largest_eigenvec[0]*np.sqrt(largest_eigenval), largest_eigenvec[1]*np.sqrt(largest_eigenval), **quiveropts)
    plt.quiver(mean_x, mean_y, smallest_eigenvec[0]*np.sqrt(smallest_eigenval), smallest_eigenvec[1]*np.sqrt(smallest_eigenval), **quiveropts)

    return


def PlotFitness(walk_folder):
    pop_file = walk_folder+"/vre.csv"

    #load data
    pop_data=pd.read_csv(pop_file)

    title_bnk="Fitness evolution of a single walk"
    destination_file=walk_folder+"/Fitness"
    ax_bnk=plt.gca()

    ax_bnk.set_title(title_bnk)
    ax_bnk.set_xlabel("Generation")
    ax_bnk.set_ylabel("Fitness")

    ax_bnk.plot(pop_data["Fitness"])

    plt.tight_layout()

    #plt.show()
    plt.savefig(destination_file+PGF_EXTENSION)
    plt.savefig(destination_file+PNG_EXTENSION)
    plt.close()

    return pop_data["Fitness"]


def PlotFitnessSample(df_fitness_data, sample_folder, max_generations):
    #destination_file=sample_folder+"/BE_plot"+PLOT_EXTENSION
    destination_file=sample_folder+"/Fitness"

    last_fitness=dict()

    best_fitness=[]
    mean_fitness=[]

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

    ax_be = plt.gca()

    ax_be.plot(mean_fitness,label="mean")
    ax_be.plot(best_fitness,label="best")
    ax_be.legend()

    # Layout and titles
    ax_be.set_title("Fitness evolution of a sample")
    ax_be.set_xlabel("Generation")
    ax_be.set_ylabel("Fitness")

    #plt.show()
    plt.savefig(destination_file+PGF_EXTENSION)
    plt.savefig(destination_file+PNG_EXTENSION)
    plt.close()
    #https://timodenk.com/blog/exporting-matplotlib-plots-to-latex/

    #prepare dataframe for experiment plot

    df_result=pd.DataFrame()
    df_result["mean"]=mean_fitness
    df_result["best"]=best_fitness

    return df_result


def PlotFitnessExperiment(df_fitness_data, experiment_folder, max_generations):
    #destination_file=sample_folder+"/BE_plot"+PLOT_EXTENSION
    destination_file=experiment_folder+"/Fitness"

    best_fitness=[]
    mean_fitness=[]

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

    ax_be = plt.gca()

    ax_be.plot(mean_fitness,label="mean")
    ax_be.plot(best_fitness,label="best")
    ax_be.legend()

    # Layout and titles
    ax_be.set_title("Fitness evolution of an experiment")
    ax_be.set_xlabel("Generation")
    ax_be.set_ylabel("Fitness")

    #plt.show()
    plt.savefig(destination_file+PGF_EXTENSION)
    plt.savefig(destination_file+PNG_EXTENSION)
    plt.close()
    #https://timodenk.com/blog/exporting-matplotlib-plots-to-latex/

    return


#def PlotFitnessLandscapeSlice(input_file,plt_title):
def PlotFitnessLandscapeSliceRobustness(walk_folder,BNK_gates,BNK_inputs):

    plt_title ="BNK("+str(BNK_gates)+","+str(BNK_inputs)+") Fitness Landscape robustness slice"

    destination_file=walk_folder+plt_title

    pop_file = walk_folder+"/vre.csv"

    x_axis="Robustness"

    #load data
    pop_data=pd.read_csv(pop_file)

    #pop_data_slice = pop_data[pop_data['Symmetry']==0.125]
    pop_data_slice = pop_data[pop_data[x_axis]==0.779362543417443]

    #X = pop_data['Robustness'].values
    X = pop_data[x_axis].values
    Y = pop_data['Fitness'].values

    plt.plot(X, Y, '.', color='black')

    plt.xlabel(x_axis)
    plt.ylabel('Fitness')

    #plt_title="BNK("+str(BNK_gates)+","+str(BNK_inputs)+") Fitness Landscape"
    plt.title(plt_title)

    plt.tight_layout()

    #plt.show()
    #plt.savefig(plt_title+PLOT_EXTENSION)
    plt.savefig(destination_file+PGF_EXTENSION)
    plt.savefig(destination_file+PNG_EXTENSION)

    plt.close()

    return


def PlotFitnessLandscapeSliceSymmetry(walk_folder,BNK_gates,BNK_inputs):

    plt_title ="BNK("+str(BNK_gates)+","+str(BNK_inputs)+") Fitness Landscape symmetry slice"

    destination_file=walk_folder+plt_title

    pop_file = walk_folder+"/vre.csv"

    x_axis="Symmetry"

    #load data
    pop_data=pd.read_csv(pop_file)

    #pop_data_slice = pop_data[pop_data['Symmetry']==0.125]
    pop_data_slice = pop_data[pop_data[x_axis]==0.779362543417443]

    #X = pop_data['Robustness'].values
    X = pop_data[x_axis].values
    Y = pop_data['Fitness'].values

    plt.plot(X, Y, '.', color='black')

    plt.xlabel(x_axis)
    plt.ylabel('Fitness')

    #plt_title="BNK("+str(BNK_gates)+","+str(BNK_inputs)+") Fitness Landscape"
    plt.title(plt_title)

    plt.tight_layout()

    #plt.show()
    #plt.savefig(plt_title+PLOT_EXTENSION)
    plt.savefig(destination_file+PGF_EXTENSION)
    plt.savefig(destination_file+PNG_EXTENSION)

    plt.close()

    return




def main(args):
    #Maximum number of plots to generate (we don't need them only, only a few to show the results
    MAX_WALK_PLOTS=50#100
    MAX_SAMPLE_PLOTS=5#25

    global PLOT_EXTENSION
    credentials=args.credentials
    experiment_id=args.experiment
    root_folder=args.folder
    use_pgf=args.latex

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

    sample_measures=dict()
    walk_measures=dict()
    sample_fitness=dict()
    walk_fitness=dict()
    experiment_type=ExperimentType.UNKNOWN
    track_variance=False

    experiment_folder=BuildExperimentFolderStr(root_folder,experiment_id)

    experiment_details=ValidateExperimentType(mysql_config,experiment_id)
    exp_type_str=experiment_details["type"].iloc[0]
    exp_track_variance_int=experiment_details["track_variance"].iloc[0]
    problem_name_str=experiment_details["name"].iloc[0]
    BNK_gates=experiment_details["bnk_gates"].iloc[0]
    BNK_inputs=experiment_details["bnk_inputs"].iloc[0]
    max_generations=experiment_details["generations"].iloc[0]


    generate_BNK_STD=False
    #Check if BNK
    match = re.search('bnk', problem_name_str, re.IGNORECASE)
    if match:
      generate_BNK_STD=True

    if exp_type_str == 'vre':
        experiment_type=ExperimentType.VRE
    elif exp_type_str == 'br':
        experiment_type=ExperimentType.BR
    elif exp_type_str == 'be':
        experiment_type=ExperimentType.BE

    if exp_track_variance_int == 1:
        track_variance=True

    df_samples=RetrieveSamples(mysql_config,experiment_id)

    #Retrieve all sample data
    sample_cnt=0
    walk_cnt=0
    for (sample_id_arr) in zip(df_samples["ID"]):
        sample_id=sample_id_arr[0]
        sample_folder=BuildSampleFolderStr(experiment_folder, sample_id)
        #if experiment_type == ExperimentType.VRE:
        #    sample_measures[sample_id]=RetrieveVRESampleData(mysql_config,sample_id)
        #elif experiment_type == ExperimentType.BR:
        #    sample_measures[sample_id]=RetrieveBRSampleData(mysql_config,sample_id)
        #elif experiment_type == ExperimentType.BE:
        #    sample_measures[sample_id]=RetrieveBESampleData(mysql_config,sample_id)


        df_walks=RetrieveWalks(mysql_config,sample_id)
        for (walk_id_arr) in zip(df_walks["ID"]):
            walk_id=walk_id_arr[0]
            walk_folder=sample_folder+"/walk_"+str(walk_id)
            
            ##Build the plots for each walk
            #if experiment_type == ExperimentType.VRE:
            #    #walk_measures[walk_id]=RetrieveVREWalkData(mysql_config,walk_id)
            #    #if walk_cnt < MAX_WALK_PLOTS: PlotVREWalk(walk_measures[walk_id], walk_folder)
            #elif experiment_type == ExperimentType.BR:
            #    #walk_measures[walk_id]=RetrieveBRWalkData(mysql_config,walk_id)
            #    #walk_fitness[walk_id]=PlotFitness(walk_folder)
            #    #if walk_cnt < MAX_WALK_PLOTS: PlotBRWalk(walk_measures[walk_id], walk_folder)
            #    #if(track_variance):
            #    #    if walk_cnt < MAX_WALK_PLOTS: PlotBRWalkVariance(walk_folder)
            #    #    if walk_cnt < MAX_WALK_PLOTS: PlotBRWalkEllipseEvolution(walk_folder, 150)
            #elif experiment_type == ExperimentType.BE:
            #    #walk_measures[walk_id]=RetrieveBEWalkData(mysql_config,walk_id)
            #    #walk_fitness[walk_id]=PlotFitness(walk_folder)
            #    #if(track_variance):
            #    #    if walk_cnt < MAX_WALK_PLOTS: 
            #    #        PlotBEWalkVariance(walk_folder)
            #    #        PlotBNK_BE(walk_folder)

            if(generate_BNK_STD):
                if walk_cnt < MAX_WALK_PLOTS:
                    #PlotBNKWalkSTD(walk_folder,BNK_gates,BNK_inputs)
                    #PlotFitnessStats(walk_folder,BNK_gates,BNK_inputs)
                    if experiment_type == ExperimentType.BE:
                        #If BNK, try to build the fitness landscape using sym X rob (check if columns in population file)
                        PlotFitnessLandscapeSliceRobustness(walk_folder,BNK_gates,BNK_inputs)
                        PlotFitnessLandscapeSliceSymmetry(walk_folder,BNK_gates,BNK_inputs)
                        PlotBNKLandscape(walk_folder,BNK_gates,BNK_inputs)
                        PlotFitnessLandscapeMultiModality(walk_folder,BNK_gates,BNK_inputs)
                        #PlotFitnessLandscapeSlice(input_file,plt_title)


            walk_cnt += 1
           
        ##Build the plots for each sample
        #if experiment_type == ExperimentType.VRE:
        #    if(track_variance):
        #        if sample_cnt < MAX_SAMPLE_PLOTS: PlotVRESampleVariancePerNeighborhood(sample_folder)
        #        if sample_cnt < MAX_SAMPLE_PLOTS: PlotVRESampleVarianceAllNeighborhoods(sample_folder)
        #elif experiment_type == ExperimentType.BR:
        #    sample_fitness[sample_id]=PlotFitnessSample(walk_fitness, sample_folder, max_generations)
        #    if sample_cnt < MAX_SAMPLE_PLOTS:
        #        PlotBRSample(walk_measures, sample_folder, track_variance)
        #elif experiment_type == ExperimentType.BE:
        #    sample_fitness[sample_id]=PlotFitnessSample(walk_fitness, sample_folder, max_generations)
        #    if sample_cnt < MAX_SAMPLE_PLOTS: 
        #        PlotBESample(walk_measures, sample_folder, track_variance)

        sample_cnt += 1


    ##Build the plots for the experiment
    #if experiment_type == ExperimentType.VRE:
    #    PlotVREExperiment(sample_measures, experiment_folder)
    #elif experiment_type == ExperimentType.BR:
    #    PlotBRExperiment(sample_measures, experiment_folder)
    #    PlotFitnessExperiment(sample_fitness, experiment_folder, max_generations)
    #    #Should be same as PlotBRSample
    #elif experiment_type == ExperimentType.BE:
    #    PlotBEExperiment(sample_measures, experiment_folder, track_variance)
    #    PlotFitnessExperiment(sample_fitness, experiment_folder, max_generations)
    
    return

if __name__ == '__main__':
    args = parse_arguments()
    if (args.credentials != None and args.experiment != None):
        try:
            main(args)
        except Exception as e:
            print("Caught exception:" + str(e))
            print(traceback.format_exc())
    else:
        print("Missing arguments!")
