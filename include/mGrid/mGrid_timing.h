/*
Copyright (c) 1998-2018 MILC Developers and Contributors

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice, this permission notice and the following
disclaimers shall be included in all copies or substantial portions of
the Software.

Neither the name of the MILC collaboration, nor the names of its
developers or contributors may be used to endorse or promote products
derived from this Software without specific prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/
/*  END LEGAL */

class HISQTimer {
/**
 * @brief Silly timer to reduce code clutter
 * @author Curtis Taylor Peterson
 */

public: 
  std::chrono::time_point<std::chrono::system_clock> epoch;
  std::chrono::time_point<std::chrono::system_clock> start;

public: GridHISQTimer() { epoch = std::chrono::system_clock::now(); }

public:
  void tic() { start = std::chrono::system_clock::now(); }
  
  void toc(const std::string& msg) {
    auto end = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << msg << ": " << elapsed << "\n";
  }

  void toctic(const std::string& msg) { toc(msg); tic(); }

  void epoch(const std::string& msg) {
    auto now = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - epoch);
    std::cout << msg << elapsed << "\n";
  }

};