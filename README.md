### AGEseqWeb Application

#### Overview
AGEseqWeb is an application designed for the identification of target editing events from the outputs of BLAT. It aligns over thousands of reads to targets and summarizes the number of target editing events supported by the number of reads per target sequence in the provided file.

#### Components
1. **UI (Streamlit-based)**:
   - Provides a user interface created using Streamlit.
   - Allows users to upload `target.txt` and FASTQ files as input.
   - Customizable parameters include Mismatch Cutoff, WT-like Report, Min Cutoff, and Indel Report.
   - Creates a folder for each request and saves the input files uploaded by the user in that folder.
   - Displays the final result on the UI, utilizing the output from the Perl API script.
   - Provides a summary table that can be downloaded locally.
  
2. **Perl API**:
   - Backend written in Perl.
   - Calls BLAT for further processing of the uploaded files.
   - Accepts parameters in the POST method.
   - Processes the files based on the parameters selected from the UI and saves the output in a `.txt` file in the same folder as created initially.

#### Server Information
- **Hosted on server**: 172.30.18.104
- **Code location**: `/data/AGEseq`
- **Accessible link**: [AGESeqWeb](http://tsailab.gene.uga.edu/AGESeqWeb/)
- **Embedded**: `/var/www/tsailabgene/AGESeqWeb/index.html`

#### Ports
- Streamlit UI port: `8088`
- Perl API port: `8087`

#### File Locations
- **Perl script**: `AGEseq.pl`
  - Located at: `/data/AGEseq`
  - Hosts the API with port and host details.
  - Requires `blat_exe` for processing data, available at: `/data/AGEseq/blat_exe`

- **Streamlit script**: `app.py`
  - Located at: `/data/AGEseq`
  - Requires `requirements.txt` for dependencies.
  - Run using: `streamlit run app.py --server.port 8088`
 
#### Systemd Configuration
- **Perl API**:
  - Service file: `/lib/systemd/system/AGEseqAPI.service`
  - Configuration:

 ```
[Unit]
Description=AGEseq API
After=multi-user.target

[Service]
WorkingDirectory=/data/AGEseq
User=tsai-apps
Type=idle
ExecStart=/usr/bin/perl /data/AGEseq/AGEseq.pl
Restart=on-failure
KillMode=process
LimitMEMLOCK=infinity
LimitNOFILE=65535
Type=simple

[Install]
WantedBy=multi-user.target

```
   -  Follow below steps to configure:
   -  ```sudo nano /lib/systemd/system/AGEseqAPI.service``` copy, paste and save above configuration.
   -  ```sudo chmod 644 /lib/systemd/system/AGEseqAPI.service```
   -  ```sudo systemctl daemon-reload```
   -  ```sudo systemctl enable AGEseqAPI.service```
   -  ```sudo systemctl start AGEseqAPI.service```
   - ```sudo systemctl status AGEseqAPI.service```

     
- **Streamlit App**:
  - Service file: `/lib/systemd/system/AGEseqAPP.service`
  - Configuration:

```
[Unit]
Description=AGEseq App
After=multi-user.target

[Service]
WorkingDirectory=/data/AGEseq
User=tsai-apps
Type=idle
ExecStart=/data/AGEseq/venv/bin/streamlit run /data/AGEseq/app.py --server.port 8088
Restart=on-failure
KillMode=process
LimitMEMLOCK=infinity
LimitNOFILE=65535
Type=simple

[Install]
WantedBy=multi-user.target

```

#### Testing
- Upload `target.txt` and files in the `read/` directory from this repository.
- Submit to view the output summary.

#### Note
- Ensure proper installation and configuration of dependencies.
- Test the application with provided test files.
- Monitor service status using `systemctl status`.
