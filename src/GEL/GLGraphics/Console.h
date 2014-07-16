/* ----------------------------------------------------------------------- *
 * This file is part of GEL, http://www.imm.dtu.dk/GEL
 * Copyright (C) the authors and DTU Informatics
 * For license and list of authors, see ../../doc/intro.pdf
 * ----------------------------------------------------------------------- */

/**
 * @file   Console.h
 * @author Anders Wang Kristensen <awk@imm.dtu.dk>
 * @date   Fri Oct 22 18:32:06 2011
 *
 * @brief  OpenGL 'Quake'-like console
 */

#ifndef __GEL_GLGRAPHICS_CONSOLE_H__
#define __GEL_GLGRAPHICS_CONSOLE_H__

#include "../GL/glew.h"

#include <functional> //make_shared, bind
#include <map> //multimap
#include <sstream> //stringstream
#include <memory> //unique_ptr
#include <vector>
#include <cassert>

namespace GLGraphics {
    
    
    

class Console
{
public:
    Console();
    ~Console() throw();

    // keyboard and display
    void display(int scaling=1);
    void keyboard(unsigned char key);
    void key_left();
    void key_right();
    void key_home();
    void key_end();
    void key_up();
    void key_down();

    //stdio-like io
    void print(const char* buffer);
    void printf(const char* format, ...);
    void newline();

    //execute command
    void execute(const char* buffer);
    void executef(const char* format, ...);

    typedef int cmd_token;

    //0-ary
    inline cmd_token reg_cmd0(const std::string& name,
                              const std::function<void ()>& f,
                              const std::string& help);

    //1-ary
    template <typename A0>
    cmd_token reg_cmd1(const std::string& name,
                       const std::function<void (const A0&)>& f,
                       const std::string& help);

    //2-ary
    template <typename A0, typename A1>
    cmd_token reg_cmd2(const std::string& name,
                       const std::function<void (const A0&,const A1&)>& f,
                       const std::string& help);

    //3-ary
    template <typename A0, typename A1, typename A2>
    cmd_token reg_cmd3(const std::string& name,
                       const std::function<void (const A0&,
                                                 const A1&,const A2&)>& f,
                       const std::string& help);

    //N-ary
    inline cmd_token reg_cmdN(const std::string& name,
                              const std::function<
                                void (const std::vector<std::string>&)>& f,
                              const std::string& help);

    //remove
    inline void unreg_cmd(cmd_token);

    //get name/help of registrered command
    const char* get_name(cmd_token) const;
    const char* get_help(cmd_token) const;

    //helper classes
    class command;

    template <typename T>
    class variable;

private:
    //make noncopyable
    Console(Console&);
    const Console& operator=(const Console&);

    static inline void from_string(const std::string& source,
                                   std::string& target)
    {
        target = source;
    }

    template <typename T>
    static inline void from_string(const std::string& source, T& target)
    {
        std::stringstream ss(source);
        if (!(ss >> target))
        {
            std::stringstream ss;
            ss << "Cannot convert "
               << source << " to '" << typeid(T).name() << "'.";
            throw std::invalid_argument(ss.str());
        }
    }

    //from string
    template <typename T>
    static T lexical_cast(const std::string& source)
    {
        T target;
        from_string(source, target);
        return target;
    }

    //to string
    template <typename S>
    static std::string to_string(const S& source)
    {
        std::stringstream ss;
        if (!(ss << source))
        {
            //converting *to* string should never fail
            throw std::invalid_argument("Cannot convert argument to string.");
        }

        return ss.str();
    }

    class command_base
    {
    public:
        inline command_base(Console& c, const std::string& help)
            : m_console(c), m_help(help) { m_id = c.m_id_counter++; }
        virtual ~command_base() {}

        virtual void execute(const std::vector<std::string>& args) const = 0;
        virtual size_t arity() const = 0;

        inline cmd_token get_id() const { return m_id; }

        inline const char* get_help() const { return m_help.c_str(); }

    protected:
        Console& m_console;
        cmd_token m_id;
        std::string m_help;
    };

    //no need to be a template
    class command0 : public command_base
    {
    public:
        typedef std::function<void ()> function_type;

        inline command0(Console& c,
                        const std::string& h, const function_type& f)
            : command_base(c,h), m_callback(f) {}

        inline void execute(const std::vector<std::string>& args) const;
        inline size_t arity() const { return 0; }

    private:
        function_type m_callback;
    };

    template <typename A0>
    class command1 : public command_base
    {
    public:
        typedef typename std::function<void (const A0&)> function_type;

        command1(Console& c, const std::string& h, const function_type& f)
            : command_base(c,h), m_callback(f) {}

        void execute(const std::vector<std::string>& args) const;
        size_t arity() const { return 1; }

    private:
        function_type m_callback;
    };

    template <typename A0, typename A1>
    class command2 : public command_base
    {
    public:
        typedef typename std::function<void (const A0&,
                                             const A1&)> function_type;

        command2(Console& c, const std::string& h, const function_type& f)
            : command_base(c,h), m_callback(f) {}

        void execute(const std::vector<std::string>& args) const;
        size_t arity() const { return 2; }

    private:
        function_type m_callback;
    };

    template <typename A0, typename A1, typename A2>
    class command3 : public command_base
    {
    public:
        typedef typename std::function<void (const A0&,
                                             const A1&,
                                             const A2&)> function_type;

        command3(Console& c, const std::string& h, const function_type& f)
            : command_base(c,h), m_callback(f) {}

        void execute(const std::vector<std::string>& args) const;
        size_t arity() const { return 3; }

    private:
        function_type m_callback;
    };

    //no need to be a template
    class commandN : public command_base
    {
    public:
        typedef std::function<
            void (const std::vector<std::string>&)> function_type;

        inline commandN(Console& c,
                        const std::string& h, const function_type& f)
            : command_base(c,h), m_callback(f) {}

        inline void execute(const std::vector<std::string>& args) const;
        inline size_t arity() const { return any_arity; }

    private:
        function_type m_callback;
    };

    //multi, so we can have multiple cmds with same name (but different arity)
    typedef std::multimap<std::string,
                          std::unique_ptr<command_base> > command_map_t;

    cmd_token add_command(const std::string& name,
                           std::unique_ptr<command_base>&& ptr);
    void remove_command(cmd_token);

    void tab_completion();

    std::vector<std::string> parse_cmdline(const char* buffer) const;

    //builtin commands
    void help();
    void help(const std::string&);
    void clear();
    void history();

    void load_history();
    void save_history() const;
    void clear_history();

    enum { any_arity = 0xFFFF };

    //draw commands
    void draw_text(int scale, int x, int y,
                   float r, float g, float b,
                   const char* buffer);

    void draw_textf(int scale, int x, int y,
                    float r, float g, float b,
                    const char* fmt, ...);
    //state
    command_map_t m_commands;
    std::vector<std::string> m_buffer;

    size_t m_history_index;
    std::vector<std::string> m_history;

    std::string::size_type m_caret;
    std::string m_current_command;

    int m_id_counter;
    bool m_is_executing;

    GLuint m_font;

    static const unsigned char g_png_data[];
    static const size_t g_png_size;
};

//0-ary
Console::cmd_token Console::reg_cmd0(const std::string& name,
                                     const std::function<void ()>& f,
                                     const std::string& help)
{
    typedef command0 command_type;
    return add_command(name,
        std::unique_ptr<command_type>(
            new command_type(std::ref(*this), help, f)));
}

void Console::command0::execute(const std::vector<std::string>& args) const
{
    assert(args.size() == arity());
    m_callback();
}

//1-ary
template <typename A0>
Console::cmd_token Console::reg_cmd1(const std::string& name,
                                     const std::function<void (const A0&)>& f,
                                     const std::string& help)
{
    typedef command1<A0> command_type;
    return add_command(name,
        std::unique_ptr<command_type>(
            new command_type(std::ref(*this), help, f)));
}

template <typename A0>
void Console::command1<A0>::execute(const std::vector<std::string>& args) const
{
    assert(args.size() == arity());
    m_callback(lexical_cast<A0>(args[0]));
}

//2-ary
template <typename A0, typename A1>
Console::cmd_token Console::reg_cmd2(const std::string& name,
                                     const std::function<void (const A0&,
                                       const A1&)>& f,
                                     const std::string& help)
{
    typedef command2<A0,A1> command_type;
    return add_command(name,
        std::unique_ptr<command_type>(
            new command_type(std::ref(*this), help, f)));
}

template <typename A0,typename A1>
void Console::command2<A0,A1>::execute(
    const std::vector<std::string>& args) const
{
    assert(args.size() == arity());
    m_callback(lexical_cast<A0>(args[0]),
               lexical_cast<A1>(args[1]));
}

//3-ary
template <typename A0, typename A1, typename A2>
Console::cmd_token Console::reg_cmd3(const std::string& name,
                                     const std::function<void (const A0&,
                                       const A1&,const A2&)>& f,
                                     const std::string& help)
{
    typedef command3<A0,A1,A2> command_type;
    return add_command(name,
        std::unique_ptr<command_type>(
            new command_type(std::ref(*this), help, f)));
}

template <typename A0,typename A1, typename A2>
void Console::command3<A0,A1,A2>::execute(
    const std::vector<std::string>& args) const
{
    assert(args.size() == arity());
    m_callback(lexical_cast<A0>(args[0]),
               lexical_cast<A1>(args[1]),
               lexical_cast<A2>(args[2]));
}

//N-ary
Console::cmd_token Console::reg_cmdN(const std::string& name,
    const std::function<void (const std::vector<std::string>&)>& f,
    const std::string& help)
{
    typedef commandN command_type;
    return add_command(name,
        std::unique_ptr<command_type>(
            new command_type(std::ref(*this), help, f)));
}

void Console::commandN::execute(const std::vector<std::string>& args) const
{
    m_callback(args);
}

void Console::unreg_cmd(cmd_token id)
{
    remove_command(id);
}

class Console::command
{
public:
    inline command() : m_console(NULL) {}

    inline void reg(Console& cs,
        const std::string& name,
        const std::function<void ()>& function,
        const std::string& help)
    {
        assert(!m_console);
        m_console = &cs;
        m_id = m_console->reg_cmd0(name, function, help);
    }

    template <typename A0>
    void reg(Console& cs,
        const std::string& name,
        const std::function<void (const A0&)>& function,
        const std::string& help)
    {
        assert(!m_console);
        m_console = &cs;
        m_id = m_console->reg_cmd1<A0>(name, function, help);
    }

    template <typename A0, typename A1>
    void reg(Console& cs,
        const std::string& name,
        const std::function<void (const A0&, const A1&)>& function,
        const std::string& help)
    {
        assert(!m_console);
        m_console = &cs;
        m_id = m_console->reg_cmd2<A0,A1>(name, function, help);
    }

    template <typename A0, typename A1, typename A2>
    void reg(Console& cs,
        const std::string& name,
        const std::function<void (const A0&,
        const A1&, const A2&)>& function,
        const std::string& help)
    {
        assert(!m_console);
        m_console = &cs;
        m_id = m_console->reg_cmd3<A0,A1,A2>(name, function, help);
    }

    inline ~command()
    {
        if (m_console)
            m_console->unreg_cmd(m_id);
    }

    inline const char* get_name() const
    {
        assert(m_console);
        return m_console->get_name(m_id);
    }

    inline const char* get_help() const
    {
        assert(m_console);
        return m_console->get_help(m_id);
    }

    inline Console* get_console() const { return m_console; }
    inline cmd_token get_id() const { assert(m_console); return m_id; }

private:
    command(command&);
    const command& operator=(const command&);

    Console* m_console;
    cmd_token m_id;
};

template <typename T>
class Console::variable
{
public:
    variable(const T& initial_value = T())
        : m_value(initial_value) {}

    void reg(Console& cs,
        const std::string& name,
        const std::string& help)
    {
        if(m_set_cmd.get_console() == 0)
        {
            m_print_cmd.reg(cs, name,
                            std::bind(&variable::print_value, this), help);
            
            m_set_cmd.reg<T>(cs, name,
                             std::bind(&variable::set_value, this, std::placeholders::_1),
                             help);
        }
    }

    const variable& operator=(const T& value) { m_value = value; return *this; }

    operator const T&() const { return m_value; }

    const char* get_name() const { return m_print_cmd.get_name(); }
    const char* get_help() const { return m_print_cmd.get_help(); }

private:
    variable(const variable&);
    const variable& operator=(const variable&);

    void print_value()
    {
        m_print_cmd.get_console()->printf("%s = %s",
            m_print_cmd.get_name(),
            Console::to_string(m_value).c_str());
    }

    void set_value(const T& value)
    {
        m_value = value;
        m_print_cmd.get_console()->execute(m_print_cmd.get_name());
    }

    T m_value;

    command m_print_cmd;
    command m_set_cmd;
};

}

#endif //__GEL_GLGRAPHICS_CONSOLE_H__
