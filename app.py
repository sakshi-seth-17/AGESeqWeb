#streamlit run app.py --server.port 8088'''
import streamlit as st
import os
import datetime
import time
import requests
import pandas as pd
from streamlit_js_eval import streamlit_js_eval
#import logging


# Configure logging
#log_format = "%(asctime)s - %(levelname)s - %(message)s"
#logging.basicConfig(filename='streamlit_logs.log', level=logging.INFO, format=log_format)

def apply_styles_to_df(value):
    return 'background-color: #f2f2f2; color: black' if value % 2 == 0 else 'background-color: white; color: black'




def save_files(files):
    # Create a folder with the current timestamp
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    folder_path = os.path.join(os.getcwd(), timestamp)
    os.makedirs(folder_path)

    # Create 'reads' subfolder for .fastq files
    reads_folder_path = os.path.join(folder_path, 'reads')
    os.makedirs(reads_folder_path)

    flag_reads = 0
    flag_targets = 0
    # Save each uploaded file
    for file in files:
        file_path = os.path.join(folder_path, file.name)
        if file.name.endswith('.fastq') or file.name.endswith('.fa'):
            file_path = os.path.join(reads_folder_path, file.name)
            flag_reads = 1
            
        if "targets.txt" in file_path:
            flag_targets = 1
            
        with open(file_path, 'wb') as f:
            f.write(file.getvalue())
            
    if (flag_targets==0 and flag_reads==0) or (flag_targets==0 and flag_reads==1) or (flag_targets==1 and flag_reads==0):
        st.error("Please upload both targets.txt and .fastq files")
        return 0
    else:
        st.success("Files uploaded successfully")
        return timestamp


    
    
# Function to send API request
def send_api_request(folder_name, mismatch_cutoff, min_cutoff, wt_like_report, indel_report):
    
    #print("----mismatch_cutoff---",mismatch_cutoff)
    #print("---min_cutoff----",min_cutoff)
    #print("---wt_like_report----",wt_like_report)
    #print("---indel_report----",indel_report)
    
    #api_endpoint = 'http://tsailab.gene.uga.edu:8502/ageseq?directory='+folder_name
    try:
        #data = requests.get('http://tsailab.gene.uga.edu:8087/ageseq?directory='+folder_name+'&mismatch_cutoff='+mismatch_cutoff+'&min_cutoff='+min_cutoff+'&wt_like_report='+wt_like_report+'&indel_report='+indel_report+'&remove_files='remove_files).json()
        data = requests.get('http://tsailab.gene.uga.edu:8087/ageseq?directory=' + folder_name + '&mismatch_cutoff=' + str(mismatch_cutoff) + '&min_cutoff=' + str(min_cutoff) + '&wt_like_report=' + str(wt_like_report) + '&indel_report=' + str(indel_report)).json()

    except:
        pass
    
def read_file_content(file_name):
    prev_size = -1  # Initial size of the file
    i = 0
    
    while True:
        # Get the current size of the file
        curr_size = os.path.getsize(file_name)
    
        # If the size stops changing, start reading the file
        if curr_size == prev_size:
            try:
                data = []
                with open(file_name, 'r') as f:
                    for line in f:
                        if "Sum:" in line and "AlignedTarget" not in line:
                            data.append(line.replace("Sum:", "").replace("\n", "").split('\t'))
    
                columns = ["INPUT", "Target", "AlignedTarget", "AlignedRead", "Total Hits", "Sub Hits",	"Indel or WT Hits",	"Indel or WT rate %",	"Pattern"]
                data_rows = data[:]
                df = pd.DataFrame(data_rows, columns=columns)
                st.write("Output")
                st.write(df)
                return df
            except Exception as e:
                #st.error(f"An error occurred while reading the file: {e}")
                break
    
        prev_size = curr_size
        time.sleep(1)  # Wait for 1 second before checking the file size again
    
    
def print_output_file(folder_path):
    file_name = "/data/AGEseq/"+folder_path+"/output.txt"
    with st.spinner("Generating output..."):
        df = read_file_content(file_name)
        
    return df
            

    

# Streamlit app
def main():
    print("here")
    st.set_page_config(page_title="AGEseq", page_icon=":microscope:", layout="wide")
    st.markdown(
        """
        <div style="background-color:#154730;padding:10px;border-radius:10px">
            <h1 style="color:white;text-align:center;">AGEseq Analysis</h1>
        </div>
        """,
        unsafe_allow_html=True
    )

    st.markdown(
        """
        <style>
        .container-3d {
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0px 0px 10px 0px rgba(0,0,0,0.3);
            background-color: #ffffff;
        }
        </style>
        """,
        unsafe_allow_html=True
    )
    

    col1, col2, col3 = st.columns([2, 1, 1])
    
    with col1:
        st.markdown('<h6><I>*Upload targets.txt and .fastq files only</I></h6>', unsafe_allow_html=True)
        # File upload
        uploaded_files = st.file_uploader(' ', accept_multiple_files=True)
        
    with col2:
        st.markdown('<h6><b>Mismatch Cutoff</b></h6>', unsafe_allow_html=True)
        mismatch_cutoff = st.slider("mismatch rate to filter low quality alignment, default 0.1 (10%)", min_value=0.0, max_value=1.0, value=0.1, step=0.01)
        st.markdown('<h6><b>Min Cutoff</b></h6>', unsafe_allow_html=True)
        min_cutoff = st.slider("cutoff to filter reads with low abundance, default 0", min_value=0, max_value=100, value=0)
        
    with col3:
        st.markdown('<h6><b>WT-like Report</b></h6>', unsafe_allow_html=True)
        wt_like_report = st.slider("report top xx WT like records, default 20", min_value=0, max_value=100, value=20)
        st.markdown('<h6><b>Indel Report</b></h6>', unsafe_allow_html=True)
        indel_report = st.slider("report top xx records with indel, default 50", min_value=0, max_value=100, value=50)
        #remove_files = st.checkbox("Remove Files", value=True)
        
    if st.button('Refresh'):
        streamlit_js_eval(js_expressions="parent.window.location.reload()")
    
    if uploaded_files:
        # Save files button
        if st.button('Process Files'):
        
            folder_name = save_files(uploaded_files)
            if folder_name != 0:                
                send_api_request(folder_name, mismatch_cutoff, min_cutoff, wt_like_report, indel_report)
                df = print_output_file(folder_name)
                st.success("Analysis completed!")
                
                
                with open(f'/data/AGEseq/{folder_name}/output.txt', "r") as file:
                    content = file.read()
                    
                #st.link_button(label="Download Output.txt", url=f'/data/AGEseq/{folder_name}/output.txt')
                st.download_button(label="Download Output.txt", data=content, file_name="Output.txt", mime="text/plain")
                st.markdown("<h6><b><I>After downloading your output.txt file, please note that the summary will no longer be visible. If you wish to view the summary again, simply click on the 'Process Files' button.</I></b></h6>", unsafe_allow_html=True)



if __name__ == '__main__':
    main()

