#--------A brief Introduction to 4mC Server---------#

The epigenetic modification, DNA N4 methylcytosine(4mC) plays an important role in DNA expression, repair, and replication. It simply plays a crucial role in restriction-modification systems. The better and accurate prediction of 4mC sites in DNA is much-needed work to understand their functional behavior that leads to help in both drug discovery and biomedical research.

For the convenience of users, particularly for the vast majority of experimental scientists, a publicly accessible web-server for 4mC-RF has been established. Moreover, to maximize users' convenience, given below a step-by-step guide on how to use it to get the desired results without the need to through the above mathematical details. The server can be accessed through http://www.biopred.org/4mC

#--------A Quick Guide to use the 4mC-RF Server---------#

Step 1.Click on the Server link from where you can either type or copy and paste the query sequences into the input box at the center of 4mC-RF Server Page. The input sequence should be in simple sequence format or in the FASTA format. A potential sequence in FASTA format consists of a single initial line beginning with the symbol ">" in the first column, followed by lines of sequence data in which nucleotides are represented using single-letter codes. Except for the mandatory symbol ">", all the other characters in the single initial line are optional and only used for the purpose of identification and description. The sequence ends if another line starting with the symbol ">" appears; this indicates the start of another sequence.

Step 2. Click on the Submit button to see the predicted results. The predicted results will be shown on the screen. The system will tell whether the DNA sequence is an 4mC positive , a non-4mC negative or an Invalid sequence. All these predicted results are fully consistent with experimental observations.

Note: You can download the Benchmark Datasets by clicking the Data link provided on the navigation bar.

############ Some remarks on used python code files ###########

The script in app.py is critical for controlling the overall functionality of the webserver. It includes all the necessary procedures that are used to interact with the user and perform operations requested upon input and navigations.

The script in extractFeatures.py includes all the important procedures that compute the features from the given DNA sequences and make predictions based on the trained models using "Model.pkl" file. Furthermore, this python script also includes all the implementations of processes used to implement statistical moments. 

passenger_wsgi.py python script file is the main application loader file for accessing the webserver.

The "templates" directory includes all the webpages in HTML that are used throughout during the application processing.

The "static" directory includes the datasets that are used in training and testing the computational model. It also includes the python packages for the current model's implementation. Most of these models are implemented and tested using Scikit-Learn library for python.

The "requirements.txt" file includes the list of all the python pips that are required to install and run the current project code.

For more information and queries kindly contact: yaser.khan@umt.edu.pk or mrfajararshad@gmail.com
