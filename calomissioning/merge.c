//
// Created by William Quinn on 18/11/2020.
//
void merge()
{
    TChain chain("T");
    for (int i=434;i<438;i++)
    {
        string file=Form("run_%d.root",i); // That should work if you have files called Sensitivity.root in subdirectories called  dir_1, dir_2 ... dir_10 Change it to match what you have!
        chain.Add(file.c_str());
    }
    chain.Merge("output.root"); // Set this to the name you want to output
}

