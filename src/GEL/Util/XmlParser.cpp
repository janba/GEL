/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

#include <fstream>
#include <string>
#include <sstream>
#include <cstring>
#include <algorithm>
#include <list>

#include "XmlParser.h"
#include "string_utils.h"

namespace Util
{
    using std::string;
    using std::map;
    using std::list;
    using std::pair;
    
    using std::ifstream;
    using std::ostream;
    using std::ostringstream;
    using std::ios_base;
    
    using std::max;
    
    using std::cerr;
    using std::endl;
    
    /////////////////////////////////////////////////////////////////
    // String handling
    /////////////////////////////////////////////////////////////////
    
    void parse_attribs(const string& s, map<string, string>& result)
    {
        string name1, name2;
        list<string> atts;
        trim_split(s, atts, "=");
        if(atts.empty())
            return;
        
        get_last(atts.front(), name1);
        list<string>::iterator i = atts.begin();
        list<string>::iterator end = atts.end();
        if(i == --end)
            return;
        
        for(++i; i != end; ++i)
        {
            get_last(*i, name2);
            result[trim(name1)] = trim(*i, " \"");
            name1 = name2;
        }
        result[trim(name1)] = trim(*end, " \"");
    }
    
    /////////////////////////////////////////////////////////////////
    // File handling
    /////////////////////////////////////////////////////////////////
    
    ifstream& seek_string(ifstream& in, const string& s, const size_t bufsize = 100)
    {
        const int bsize = static_cast<int>(max(s.size(), bufsize));
        const int n = static_cast<int>(s.size());
        char* buf = new char[bsize + 1];
        char s0 = s[0];
        
        in.get(buf, bsize, s0);
        in.clear();
        in.read(buf, n);
        buf[n] = '\0';
        while(in && strcmp(s.c_str(), buf) != 0)
        {
            in.get(buf, bsize, s0);
            in.clear();
            in.read(buf, n);
            buf[n] = '\0';
        }
        
        delete [] buf;
        return in;
    }
    
    ifstream& read_until(ifstream& in, string& s_in, const string s, const size_t bufsize = 500)
    {
        const int bsize = static_cast<int>(max(s.size(), bufsize));
        const int n = static_cast<int>(s.size());
        char* buf = new char[bsize + 1];
        char s0 = s[0];
        ostringstream ostr;
        
        in.get(buf, bsize, s0);
        ostr << buf;
        in.clear();
        in.read(buf, n);
        buf[n] = '\0';
        while(in && strcmp(s.c_str(), buf) != 0)
        {
            ostr << buf;
            in.get(buf, bsize, s0);
            ostr << buf;
            in.clear();
            in.read(buf, n);
            buf[n] = '\0';
        }
        s_in = ostr.str();
        
        delete [] buf;
        return in;
    }
    
    ifstream& operator>>(ifstream& in, XmlHead& fhead)
    {
        seek_string(in, "<?xml");
        
        string head;
        read_until(in, head, "?>");
        
        fhead.is_xml = in.good();
        parse_attribs(head, fhead.atts);
        
        return in;
    }
    
    ifstream& operator>>(ifstream& in, XmlElement& elem)
    {
        seek_string(in, "<");
        
        string head;
        read_until(in, head, ">");
        if(head.empty() || head[0] == '!') return in;
        if(head[0] == '/')
        {
            if(elem.parent && elem.parent->name == head.substr(1))
                in.setstate(ios_base::eofbit);
            else
                elem.name = "";
            return in;
        }
        
        bool has_body = true;
        
        if(head[head.size() - 1] == '/')
        {
            has_body = false;
            head.erase(head.size() - 1);
        }
        
        get_first(head, elem.name);
        elem.name = trim(elem.name);
        parse_attribs(head, elem.atts);
        
        if(has_body)
        {
            delete elem.body;
            elem.body = new XmlBody(elem.doc);
            elem.body->element.parent = &elem;
            read_until(in, elem.body->text, "<");
            in.putback('<');
        }
        return in;
    }
    
    ifstream& operator>>(ifstream& in, XmlBody& body)
    {
        read_until(in, body.text, "<");
        in.putback('<');
        return in;
    }
    
    /////////////////////////////////////////////////////////////////
    // Methods
    /////////////////////////////////////////////////////////////////
    
    XmlElement::~XmlElement()
    {
        delete body;
    }
    
    void XmlBody::process_element()
    {
        if(!doc) return;
        if((doc->infile >> element) && !doc->infile.eof())
        {
            XmlElementHandler h = doc->handlers[element.name];
            if(h)
                h(element);
            else
                element.process_elements();
        }
    }
    
    void XmlElement::process_elements()
    {
        if(!doc) return;
        if(!body) return;
        
        while((doc->infile >> body->element) && !doc->infile.eof())
        {
            XmlElementHandler h = doc->handlers[body->element.name];
            if(h)
                h(body->element);
            else
                body->element.process_elements();
        }
        doc->infile.clear();
    }
    
    XmlDoc::XmlDoc(const char *filename)
    : infile(filename), body(0)
    {
        if(!infile)
        {
            cerr << "cannot open input file" << filename << endl;
            return;
        }
        
        if(infile >> head)
            infile >> body;
        else
            cerr << filename << " is not a valid xml-file" << endl;
        
        body.doc = this;
        body.element.doc = this;
    }
    
    void XmlDoc::process_elements()
    {
        while(!infile.eof() && (infile >> body.element))
        {
            XmlElementHandler h = handlers[body.element.name];
            if(h)
                h(body.element);
            else
                body.element.process_elements();
        }
        infile.clear();
    }
    
    /////////////////////////////////////////////////////////////////
    // Output handling
    /////////////////////////////////////////////////////////////////
    
    ostream& operator<<(ostream& out, const pair<string, string>& attrib)
    {
        return out << attrib.first << "=\"" << attrib.second << "\"";
    }
    
    ostream& operator<<(ostream& out, const XmlHead& head)
    {
        out << "<?xml";
        for(map<string, string>::const_iterator i = head.atts.begin(); i != head.atts.end(); ++i)
            out << " " << *i;
        out << "?>";
        return out;
    }
}