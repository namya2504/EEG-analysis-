
How to Download the Data:

1. Visit the PhysioNet Sleep-EDF Database:
   https://physionet.org/content/sleep-edfx/1.0.0/

2. Download the PSG files (polysomnography recordings):
   - Navigate to `sleep-cassette/`
   - Download any file: `SC4001E0-PSG.edf`, `SC4002E0-PSG.edf`, etc.

3. Place the `.edf` files in the folder named `data/`

Sleep EEG Analysis Model

This project represents my first hands-on experience with computational neuroscience and signal processing. I taught myself MATLAB to analyze sleep EEG data and explore how different brain wave frequencies change throughout the night. 
I focused on understanding basic signal processing concepts and frequency analysis rather than building a predictive model.


Explanation of the structure:

Step 1. Data Loading & Exploration
- Learned the EDF file format and wrote a basic reader
- Loaded sleep polysomnography recordings
- Created initial visualizations to understand the data

Step 2. Preprocessing
- Applied bandpass filter (0.5-32 Hz) to remove noise
- Segmented continuous recording into 10-second epochs
- Prepared clean data for analysis

Step 3. Frequency Analysis
- Extracted power in each frequency band using Welch's method
- Calculated band power for all segments
- Visualized how frequency patterns change over time

Step 4. Pattern Exploration
- Compared first half vs second half of night
- Identified unusual segments 
- Explored relationships between frequency bands




Skills acquired while working on this project: 
- MATLAB programming basics
- Signal processing (filtering, spectral analysis)
- Frequency analysis
- Data visualization 


Skill I need to learn:
- Time-frequency analysis (wavelets) for better temporal resolution
- Advanced feature extraction (entropy measures, complexity metrics)
- Machine learning for classification
- Statistical hypothesis testing and validation
- Working with multi-subject datasets

## Limitations

Data:
- Single subject
- Memory constraints
- No sleep stage labels 
- Healthy controls only

Methods:
- Basic power analysis
- Simplified approach
- No statistical tests



## Repository Structure
