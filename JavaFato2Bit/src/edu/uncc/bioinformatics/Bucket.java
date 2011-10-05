package edu.uncc.bioinformatics;

public class Bucket {
 
	public Bucket( boolean f){
		new_seq =f;
	}
	
	public void setDone( boolean set){
		done = set;
	}
	public boolean getDone(){
		return done;
	}
	
    public Bucket( char fa[] , String fa_name){
        this.fa = fa;
        faName = fa_name;
        new_seq = false;
        done = false;
    }
    
    char fa[];
    boolean new_seq;
    boolean done;
    String faName;

}
