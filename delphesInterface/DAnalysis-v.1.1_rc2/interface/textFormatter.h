#ifndef TEXTFORMATTER_H_
#define TEXTFORMATTER_H_
/** \class d_ana::textFormatter
 *
 * class to format text
 *
 * \original author: Jan Kieseler
 *
 * more docu
 *
 */


#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include "TString.h"

namespace d_ana{

class textFormatter {

public:

	textFormatter() :trim_(" \t"),comment_("#"),delimiter_(","){}

	/**
	 * sets all chars that are used for trimming
	 * can be more than one, e.g. "\t\n "
	 */
	void setTrim(const std::string& tr){trim_=tr;}
	/**
	 * defines one char that will serve as comment indicator.
	 * all text following that char will be ignored until the next line starts
	 */
	void setComment(const std::string& c){comment_=c;}
	/**
	 * sets a delimiter for individual entries per line (e.g. "," or " ")
	 */
	void setDelimiter(const std::string& d){delimiter_=d;}

	std::string & ltrim(std::string & str) const;
	std::string & rtrim(std::string & str) const;
	std::string & trim(std::string & str) const;
	std::string & trimcomments(std::string & str) const;


	/**
	 * returns formatted text fragments according to the farmet rules defined previously
	 */
	template<class T>
	std::vector<T> getFormatted(const std::string& in)const{
		std::vector<T> out;
		using namespace std;
		string s=in;
		trimcomments(s);
		trim(s);
		istringstream ss(s);
		while (ss)
		{
			string s2;
			if (!getline( ss, s2, *delimiter_.data() )) break;
			if(debug)
				std::cout << "got \"" << s2 << "\""<<std::endl;
			trim(s2);
			if(debug)
				std::cout << "trimmed to \"" << s2 << "\"" << std::endl;
			if(s2.length()>0)
				out.push_back( (T)s2 );
		}
		return out;
	}

	std::vector<std::string> getFormatted(const std::string& in)const{
		return getFormatted<std::string>(in);
	}
	/**
	 * returns file name with extension but without path
	 */
	static std::string getFilename(const std::string& pathtofile, bool withoutextension=false);
	/**
	 * returns file extension or empty string if none
	 */
	static std::string getFileExtension(const std::string& pathtofile);
	/**
	 * strips file of extension
	 */
	static std::string stripFileExtension(const std::string& pathtofile);
	/**
	 * strips file of dir
	 */
	static std::string stripFileDir(const std::string& pathtofile);
	/**
	 * gets file directory
	 */
	static std::string getFileDir(const std::string& pathtofile);

	/**
	 * adds a suffix to file name but keeps extension, e.g. path/bla.txt is transformed to
	 * path/bla_<suffix>.txt
	 */
	static std::string addFilenameSuffix(const std::string& pathtofile, const std::string& suffix);

	/**
	 * left offset does not count to the number of max chars before linebreak
	 */
	static std::string makeCompatibleFileName(const std::string &);

	static std::string splitIntoLines(const std::string &,const size_t& maxchars, const size_t& leftoffset, size_t skipoffset=0);

	template<class T>
	static std::string toString(T in) {
		std::ostringstream s;
		s << in;
		std::string out = s.str();
		return out;
	}

	static std::string fixLength(const char* in, size_t l, bool truncate=true);
	static std::string fixLength(const TString & in, size_t l, bool truncate=true);
	/**
	 * returns a string of fixed length
	 */
	static std::string fixLength(const std::string &, size_t, bool truncate=true);

	/**
	 * switch for more output
	 */
	static bool debug;

protected:
	std::string trim_,comment_,delimiter_;


};

}




#endif /* TEXTFORMATTER_H_ */
