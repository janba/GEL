/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file XmlParser.h
 * @brief Simple XML parser.
 */
#ifndef __UTIL_XMLPARSER_H
#define __UTIL_XMLPARSER_H

#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <map>

namespace Util
{
  struct XmlHead
  {
    XmlHead() : is_xml(false) { }

    bool is_xml;
    std::map<std::string, std::string> atts;
  };
  
  struct XmlBody;
  class XmlDoc;

  struct XmlElement
  {
    XmlElement(XmlDoc* in_doc, XmlElement* parent_elem) : body(0), doc(in_doc), parent(parent_elem) { }
    XmlElement(XmlDoc* in_doc) : body(0), doc(in_doc), parent(0) { }
    ~XmlElement();

    void process_elements();

    std::string name;
    std::map<std::string, std::string> atts;
    XmlBody* body;
    XmlDoc* doc;
    XmlElement* parent;
  };

  struct XmlBody
  {
    XmlBody(XmlDoc* parent) : element(parent), doc(parent) { }

    void process_element();

    std::string text;
    XmlElement element;
    XmlDoc* doc;
  };

  typedef void (*XmlElementHandler)(XmlElement&);

  class XmlDoc
  {
  public:
    XmlDoc(const char* filename);

    bool is_valid() const { return head.is_xml; }
    void add_handler(const std::string& element_name, const XmlElementHandler& h) { handlers[element_name] = h; }
    void process_elements();
    void close() { infile.close(); }

    XmlHead head;

  private:
    friend struct XmlElement;
    friend struct XmlBody;

    std::ifstream infile;
    XmlBody body;
    std::map<std::string, XmlElementHandler> handlers;
  };

  std::ostream& operator<<(std::ostream& out, const XmlHead& head);
}

#endif // XMLPARSER_H