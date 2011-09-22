package edu.uncc.bioinformatics;
/** 
Copyright (C) 2011  Jeremy Villalobos

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

This program was a conversion from C languge.  The original code is from 
http://code.google.com/p/zinba/
**/
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
