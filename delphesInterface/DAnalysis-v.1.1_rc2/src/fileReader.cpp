/*
 * fileReader.cc
 *
 *  Created on: Oct 10, 2013
 *      Author: kiesej
 */
#include "../interface/fileReader.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdio>
#include <algorithm>

namespace d_ana{



bool fileReader::debug=false;


/**
 * ignores any white spaces!!
 */
void fileReader::setStartMarker( std::string d){
	d.erase(std::remove_if(d.begin(), d.end(), ::isspace), d.end());
	start_=d;
}
/**
 * ignores any white spaces!!
 */
void fileReader::setEndMarker( std::string d){
	d.erase(std::remove_if(d.begin(), d.end(), ::isspace), d.end());
	end_=d;
}




void fileReader::readFile(const std::string &filename){

	using namespace std;
	lines_.clear();
	ifstream myfile (filename.data(), ios::in);
	if (myfile.is_open())
	{
		tempinfilename_=filename;
		if(debug){
			std::cout << "fileReader::readFile: opened file.\n searching for start marker " << start_
					<<" until end marker " << end_ << "\n"
					<< "delimiter: " << delimiter_ << " comments: " << comment_ << std::endl;

		}
		if(blindmode_){
			vector <string> record;
			lines_.push_back(record);
		}
		bool started=false;
		bool ended=false;
		while ( myfile.good() )
		{
			if(ended)
				break;
			string s;

			if (!getline( myfile, s )) break;

			if(blindmode_){
				lines_.at(0).push_back(s);
				continue;
			}
			else{
				trimcomments(s);
				trim(s);
			}

			if(s.size()<1)
				continue;
			string d=s;
			d.erase(std::remove_if(d.begin(), d.end(), ::isspace), d.end());
			if(!blindmode_ && !started && start_.size() !=0 && d.find(start_) == string::npos) //wait for start marker on line basis
				continue;
			else if(!started){
				started=true;
				if(start_.size() >0){
					if(debug)
						std::cout << d << " is start marker, begin reading" << std::endl;
					continue; //dont read start marker line itself
				}
			}
			if(started && end_.size() >0 && d.find(end_) != string::npos){
				if(debug)
					std::cout << d << " is end marker, stop reading" << std::endl;
				ended=true;
				break;
			}

			if(debug && !started)
				std::cout << "found start marker in line" << std::endl;

			//if(end_.size() !=0 && s.compare(end_) == 0) //wait for start marker
			//break;

			istringstream ss( s );
			vector <string> record;
			bool noentry=true;
			while (ss)
			{

				string s2;
				if (!getline( ss, s2, *delimiter_.data() )) break;
				/*	if(!started && start_.size() !=0 && s2 != start_){ //wait for start marker on word by word basis
					if(debug)
						std::cout << s2 << " not yet startmarker" << std::endl;
					continue;
				}
					else if(!started){
					started=true;
					if(start_.size() >0){
						if(debug)
							std::cout << s2 << " is start marker, begin reading" << std::endl;
						continue; //dont read start marker line itself
					}
				}
				if(end_.size() >0 && s2 == end_){
					if(debug)
						std::cout << s2 << " is end marker, stop reading" << std::endl;
					ended=true;
					break;
				} */
				if(debug)
					std::cout << "read \"" << s2 << "\""<<std::endl;
				trim(s2);
				if(debug)
					std::cout << "trimmed to \"" << s2 << "\"" << std::endl;
				if(s2.size()<1)
					continue;
				noentry=false;
				record.push_back( s2 );
			}
			if(!ended && !noentry){
				lines_.push_back(record);
				record.clear();
			}
		}

		myfile.close();
	}
	else{
		cout << "fileReader::readFile: could not read file "<< filename <<std::endl;
		return;
	}
}

std::string fileReader::getReJoinedLine(size_t line) const{
	if(line>=lines_.size()){
		throw std::out_of_range("fileReader::getReJoinedLine");
	}
	std::string out;
	for(size_t i=0;i<lines_.at(line).size();i++){
		out+=lines_.at(line).at(i);
		if(i<lines_.at(line).size()-1)
			out+= delimiter_;
	}
	return out;
}


std::vector<std::string> fileReader::getMarkerValues(const std::string& markername)const{

	std::vector<std::string> out;
	for(size_t i=0;i<nLines();i++){
		std::string fullmarker;
		std::vector<std::string> formattedentries;
		if(delimiter_ == " " || delimiter_ == "-" || delimiter_ == "]" ||delimiter_ == "["  ){
			fullmarker=getReJoinedLine(i);
			textFormatter tf;
			fullmarker=getReJoinedLine(i);
			tf.setComment(comment_);
			tf.setDelimiter(",");
			tf.setTrim(trim_);
			formattedentries=tf.getFormatted(fullmarker);
		}
		else{
			formattedentries=lines_.at(i);
		}
		for(size_t entr=0;entr<formattedentries.size();entr++){

			textFormatter tfentr;
			tfentr.setDelimiter("-");
			tfentr.setComment(comment_);
			tfentr.setTrim(trim_);
			std::vector<std::string> varsandvals=tfentr.getFormatted(formattedentries.at(entr));

			if(varsandvals.size()<1)
				continue;
			if(varsandvals.at(0).find("[") ==  std::string::npos) //no marker
				continue;
			tfentr.setTrim("[");
			tfentr.ltrim(varsandvals.at(0));
			tfentr.setTrim(" ");
			tfentr.trim(varsandvals.at(0));
			if(varsandvals.at(0) == markername){

				tfentr.setTrim("[]");
				tfentr.trim(formattedentries.at(entr));
				tfentr.setTrim(markername);
				tfentr.trim(formattedentries.at(entr));
				tfentr.setTrim("- ");
				tfentr.trim(formattedentries.at(entr));

				out.push_back( formattedentries.at(entr));
			}

		}

	}

	return out;
}



const std::vector<std::string>& fileReader::getData(const size_t &line) const{
	if(line>=lines_.size()){
		throw std::out_of_range( "fileReader::getData: line out of range");
	}
	return lines_.at(line);
}
/**
 * if file has entries like somevariable=blabla
 * getValue("somevariable") will return "blabla"
 * if value is definied several times, an exception is thrown
 * less performance but safer
 * can be switched off by bool
 * then last value is returned
 *
 * if value is not found empty string will be returned
 *
 *
 * fixed format: commas as separators (if any)
 */
std::string fileReader::getValueString(const std::string & val, bool checkdoubles){
	std::string out="";
	size_t count=0;
	//search all entries:
	for(size_t line=0;line<nLines();line++){
		std::string thisline;
		textFormatter tf;
		std::vector<std::string> formattedentries;
		if(delimiter_!=","){
			thisline=getReJoinedLine(line);
			tf.setComment(comment_);
			tf.setDelimiter(",");
			tf.setTrim(trim_);
			formattedentries=tf.getFormatted(thisline);
		}
		else{
			formattedentries=lines_.at(line);
		}
		//split according to format, allow for commas as separators and spaces in names

		for(size_t entr=0;entr<formattedentries.size();entr++){
			textFormatter tfentr;
			tfentr.setDelimiter("=");
			tfentr.setComment(comment_);
			tfentr.setTrim(trim_);
			std::vector<std::string> varsandvals=tfentr.getFormatted(formattedentries.at(entr));
			if(varsandvals.size() < 2) continue;


			if(varsandvals.at(0) == val){ //found
				if(checkdoubles && count > 0){//except
					std::cout << "fileReader::getValue: value " << val<<  " defined twice! Source of errors!" <<std::endl;
					throw std::logic_error("fileReader::getValue: value defined twice! Source of errors!");
				}

				out=varsandvals.at(1);
				count++;
			}
		}
	}

	return out;
}



std::string fileReader::dumpFormattedToTmp()const{
	char buffer [256]="/tmp/tmpfileXXXXXX";
	mkstemp (buffer);
	std::string filename(buffer);
	std::ofstream ofs (filename, std::ofstream::out);
	for(size_t i=0;i<nLines();i++){
		for(size_t entr=0;entr<lines_.at(i).size();entr++){
			ofs << lines_.at(i).at(entr) << delimiter_;
		}
		ofs << std::endl;
	}

	ofs.close();
	return filename;
}

std::vector<std::string> fileReader::readList(const std::string& infile,
		const std::string& startmarker, const std::string& endmarker, const std::string& comment, const std::string& delimiter)const{
	fileReader fr;
	fr.setStartMarker(startmarker);
	fr.setEndMarker(endmarker);
	if(comment.length()<1)
		fr.setComment(comment_);
	else
		fr.setComment(comment);
	if(delimiter.length())
		fr.setDelimiter(delimiter);
	else
		fr.setDelimiter(delimiter_);
	fr.readFile(infile);
	size_t nlines=fr.nLines();
	std::vector<std::string> out;
	for(size_t i=0;i<nlines;i++){
		out.insert(out.end(),   fr.getData(i).begin(), fr.getData(i).end());
	}
	return out;
}

}//d_ana
