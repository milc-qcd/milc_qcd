#include <stdio.h>

#include <String.h>

int main ( int argc, char* argv [] )
{
  String* a = new_String ( 64, "This is string A" );
  String* b = new_String (  8, "This is string B which is bigger than A!" );
  String* c = StringClone ( a );

  printf ( "String* a = new_String ( 64, \"This is string A\" );\n" );
  printf ( "String* b = new_String (  8, \"This is string B which is bigger than A!\" );\n" );
  printf ( "String* c = StringClone ( a );\n" );
  dumpString ( stdout, "main::a", a );
  dumpString ( stdout, "main::b", b );
  dumpString ( stdout, "main::c", c );

  printf ( "\nStringCopy ( c, b );\n" );
  StringCopy ( c, b );
  dumpString ( stdout, "main::c", c );

  printf ( "\nStringConcat ( c, a );\n" );
  StringConcat ( c, a );
  dumpString ( stdout, "main::c", c );

  printf ( "StringConcatChars ( c, \" + something extra!\" );\n" );
  StringConcatChars ( c, " + something extra!" );
  dumpString ( stdout, "main::c", c );

  printf ( "StringCopyChars ( c, \"c after StringCopyChars.\" );\n" );
  StringCopyChars ( c, "c after StringCopyChar." );
  dumpString ( stdout, "main::c", c );

  {
    char* fmt = 
      "/ChiSq_Nx1/fit/params/state[%d]/overlap/Zed[%d]/prior/width/text()";

    printfToString ( a, fmt, 1, 3 );
    dumpString ( stdout, "main::a", a );
  }

  {
    String* d = new_String ( 12, "" );
    char* fmt = 
      "/ChiSq_Nx1/fit/params/state[%d]/overlap/Zed[%d]/prior/width/text()";

    printf ( "String* d = new_String ( 12, \"\" );\n" );
    dumpString ( stdout, "main::d", d );

    printf ( "printfToString ( d, fmt, 1, 3 );\n" );
    printfToString ( d, fmt, 1, 3 );
    dumpString ( stdout, "main::d", d );

    delete_String ( d );
  }


  delete_String ( c );
  delete_String ( b );
  delete_String ( a );

  return 0;
}
