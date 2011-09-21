package edu.uncc.bioinformatics;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

public class LineFile {
	File mFile;
	BufferedReader reader;
	String ReuseLine = null;
	public LineFile(String file_name) throws FileNotFoundException{
		mFile = new File( file_name );
		DataInputStream stream = new DataInputStream( new FileInputStream(mFile) );
		reader = new BufferedReader( new InputStreamReader( stream ) );
	}
	public String lineFileNext() throws IOException{
		String ans = null;
		if( ReuseLine == null){
			ans = reader.readLine();
		}else{
			ans = ReuseLine;
			ReuseLine = null;
		}
		return ans;
	}
	public File getFile(){
		return mFile;
	}
	public void lineFileReuse(String line){
		ReuseLine = line;
	}
	public void close() throws IOException{
		reader.close();
	}
}
