/*
 * textFormatter.cc
 *
 *  Created on: May 6, 2014
 *      Author: kiesej
 */


#include "../interface/textFormatter.h"
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <algorithm>

namespace d_ana{

bool textFormatter::debug=false;


/**
 * cuts on right side of the string
 * example:
 *     setTrim(#)
 *     string="#######this is a string####"
 *     rtrim(string)
 *     now: string="#######this is a string"
 */
std::string & textFormatter::rtrim(std::string & str) const{
	if(trim_.length()<1) return str;
	size_t endpos = str.find_last_not_of(trim_);
	if( std::string::npos != endpos )
	{
		str = str.substr( 0, endpos+1 );
	}
	else
		str.clear(); //only contains trim chars
	return str;
}
/**
 * cuts on right side of the string
 * example:
 *     setTrim(#)
 *     string="#######this is a string####"
 *     ltrim(string)
 *     now: string="this is a string####"
 */
std::string & textFormatter::ltrim(std::string & str) const{
	if(trim_.length()<1) return str;

	size_t startpos = str.find_first_not_of(trim_);
	if( std::string::npos != startpos )
	{
		str = str.substr(startpos );
	}
	else
		str.clear(); //only contains trim chars
	return str;
}
/**
 * cuts everything that follows the comment marker (one char)
 * example:
 *     setComment("&")
 *     string="this is a string &with a comment"
 *     trimcomments(string)
 *     now: string=="this is a string "
 */
std::string & textFormatter::trimcomments(std::string & str) const{
	if(comment_.length()<1) return str;
	if(str.length() <2 && str==comment_) str="";
	size_t endpos = str.find(comment_);
	if( std::string::npos != endpos)
	{
		str = str.substr( 0, endpos);
	}

	return str;
}

// trim from both ends
std::string &textFormatter::trim(std::string &s) const {
	ltrim(rtrim(s));

	return s;
}



//static member functions

std::string textFormatter::getFilename(const std::string& pathtofile, bool withoutextension){
	std::string out;
	using namespace std;
	string s=pathtofile;
	istringstream ss(s);
	while (ss)
	{
		string s2;
		if (!getline( ss, s2, *"/" )) break;
		if(debug)
			std::cout << "got \"" << s2 << "\""<<std::endl;
		if(s2.size()>0)
			out = s2;
	} //just keep last

	if(withoutextension){
		out = out.substr(0, out.find_last_of("."));
	}
	return out;
}


std::string textFormatter::getFileExtension(const std::string& pathtofile){
	std::string out=getFilename(pathtofile);
	using namespace std;
	string s=out;
	istringstream ss(s);
	while (ss)
	{
		string s2;
		if (!getline( ss, s2, *"." )) break;
		if(debug)
			std::cout << "got \"" << s2 << "\""<<std::endl;
		if(s2.size()>0)
			out = s2;
	} //just keep last

	//in case no extension
	if(out==s)
		return "";
	return out;
}
std::string textFormatter::stripFileExtension(const std::string& pathtofile){
	std::string temp=getFilename(pathtofile);
	using namespace std;
	string s=temp;
	string out="";
	istringstream ss(s);
	size_t endpos = pathtofile.find_last_of(".");
	if( std::string::npos != endpos )
	{
		out = pathtofile.substr( 0, endpos);
	}


	return out;
}
std::string textFormatter::stripFileDir(const std::string& pathtofile){
	std::string temp=getFilename(pathtofile);
	using namespace std;
	string s=temp;
	string out="";
	istringstream ss(s);
	size_t endpos = pathtofile.find_last_of("/");
	if( std::string::npos != endpos )
	{
		out = pathtofile.substr( endpos+1, pathtofile.length());
	}

	return out;
}

std::string textFormatter::getFileDir(const std::string& pathtofile){
	using namespace std;
	string str=pathtofile;
	size_t endpos = str.find_last_of("/");
	if( std::string::npos != endpos )
	{
		str = str.substr( 0, endpos+1 );
	}
	return str;

}
std::string textFormatter::addFilenameSuffix(const std::string& pathtofile, const std::string& suffix){

	std::string onlyname=stripFileExtension(pathtofile);

	std::string extension=getFileExtension(pathtofile);
	return onlyname+suffix+"."+extension;

}

std::string textFormatter::makeCompatibleFileName(const std::string &in){
	std::string out=in;
	std::replace( out.begin(), out.end(), '#', '_');
	std::replace( out.begin(), out.end(), '/', '_');
	std::replace( out.begin(), out.end(), '{', '_');
	std::replace( out.begin(), out.end(), '}', '_');
	std::replace( out.begin(), out.end(), ' ', '_');
	std::replace( out.begin(), out.end(), '\\', '_');
	std::replace( out.begin(), out.end(), '-', '_');
	return out;
}



std::string textFormatter::splitIntoLines(const std::string  &in,const size_t& maxchars, const size_t& leftoffset, size_t skipoffset){
	std::string out;
	const size_t& fullsize=in.length();

	size_t pos0=in.find_first_not_of(' ')+1,
			pos1=pos0,
			lastbreak=0;
	std::string offsetstring(leftoffset,' ');
	size_t iter=0;
	size_t nsplitchars=0;
	bool wordsplit=false;
	while(1){
		iter++;
		if(iter>=in.length()){
			//failure
			break;
		}
		if(fullsize-lastbreak < maxchars){
			if(lastbreak){
				out+="\n";
				out+=offsetstring;
				out+=in.substr(lastbreak+1,fullsize);
			}
			else{
				if(!skipoffset)
					out+=offsetstring;
				out+=in;
			}
			break;
		}

		pos0=in.find(' ',pos0+1);
		pos1=in.find(' ',pos0+1);


		if(pos0 - lastbreak < maxchars && pos1 - lastbreak >= maxchars){//add break
			nsplitchars=pos0-lastbreak;
		}
		else if(pos0-lastbreak >= maxchars-1){
			nsplitchars=maxchars-1;
			wordsplit=true;
		}

		if(nsplitchars){
			if(lastbreak)
				out+="\n";
			if(!skipoffset)
				out+=offsetstring;
			else
				skipoffset--;
			if(lastbreak)
				out+=in.substr(lastbreak+1,nsplitchars);
			else
				out+=in.substr(lastbreak,nsplitchars);

			if(wordsplit)
				out+="-";
			pos0=lastbreak+nsplitchars;
			lastbreak=pos0;
			nsplitchars=0;
			wordsplit=false;
		}


		//pos0=pos1;
	}
	return out;
}

std::string textFormatter::fixLength(const char*  in, size_t l, bool truncate){
	return fixLength(std::string(in),l,truncate);
}

std::string textFormatter::fixLength(const TString & in, size_t l, bool truncate){
	return fixLength((std::string)in.Data(),l,truncate);
}

std::string textFormatter::fixLength(const std::string & in, size_t l, bool truncate){
	if(in.length()==l)
		return in;
	if(truncate){
		if(in.length()>l)
			return std::string(in.begin(),in.begin()+l);
		else
			return in+std::string(l-in.length(),' ');
	}
	else{
		if(in.length()>l)
			return std::string(in.begin()+(in.length()-l),in.end());
		else
			return std::string(l-in.length(),' ')+in;
	}
}

}//ns
