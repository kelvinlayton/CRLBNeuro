% BibTeX bibliography style `subset'
% by David Kotz dfk@cs.dartmouth.edu
% August 1992, May 1993, February 1999
%   This takes an aux file and produces an unusual .bbl file,
% that is actually in .bib format -- but it has all of the references,
% and strings expanded. Thus, it becomes a complete, portable .bib
% file. Crossrefs are retained but redundant.
% This one doesn't include 'comment' in output.
%  modified from
   % BibTeX bibliography style `abstract' 
   % by David Kotz dfk@cs.duke.edu
   % March 1989
   %  modified from
   %  BibTeX standard bibliography style `alpha'
	 % version 0.99a for BibTeX versions 0.99a or later, LaTeX version 2.09.
      % Copyright (C) 1985, all rights reserved.
      % Copying of this file is authorized only if either
      % (1) you make absolutely no changes to your copy, including name, or
      % (2) if you do make changes, you name it something other than
      % btxbst.doc, plain.bst, unsrt.bst, alpha.bst, and abbrv.bst.
      % This restriction helps ensure that all standard styles are identical.
      % The file btxbst.doc has the documentation for this style.

% DFK added abstract, comment, keyword to standard set
% 5/93 DFK added private, ignored when subsetting
% 7/94 DFK added URL
ENTRY
  { abstract
    address
    author
    booktitle    
    category
    chapter
    comment
    earlier
    edition
    editor
    howpublished
    institution
    journal
    key
    keyword
    later
    month
    note
    number
    organization
    pages
    private
    publisher
    school
    series
    title
    type
    URL
    volume
    year
  }
  {}
  { label extra.label sort.label }

STRINGS { s t }

% output a field; the top is the text, the next to top is the name of
% the field, and we know that the text is nonnull.
% we start by finishing the previous line
FUNCTION {output.field}
{
  ", " write$		% write comma and newline
  newline$
  "  " write$		% and a respectable indentation

  swap$			% now field name is on top
  " = {" *		% put = {
  swap$ *		% now text is back on top; concatenate
  "}" *			% put }
  write$
}

% output if top not empty; below the top is the name of the field
FUNCTION {output}
{ duplicate$ empty$
    { pop$ pop$ }
    'output.field
  if$
}

% prints each field 
FUNCTION {output.bibfields}
{
% The fields, in what seems to be a good order
  "key" key output
  "author" author output
  "title" title output

  "journal" journal output
  "booktitle" booktitle output

  "chapter" chapter output
  "edition" edition output
  "editor" editor output

   crossref missing$
     { }
     { "crossref" crossref output }
   if$

  "year" year output
  "month" month output
  "series" series output
  "volume" volume output
  "number" number output
  "type" type output

  "pages" pages output

  "institution" institution output
  "school" school output
  "organization" organization output
  "publisher" publisher output
  "howpublished" howpublished output
  "address" address output

  "category" category output

  "note" note output

  "earlier" earlier output
  "later" later output

  "URL" URL output
  "keyword" keyword output
  "abstract" abstract output
%  "comment" comment output
}

% Takes name of type as argument, and prints @ line
% then it prints each field 
FUNCTION {output.bibitem}
{
% write the @ line
  "@" swap$ * write$    		  % @type
  "{" cite$ * write$  		  % {citekey

  output.bibfields

% finish entry
  newline$
  "}" write$ 
  newline$
  newline$
}

FUNCTION {article}{ "Article" output.bibitem }

FUNCTION {book}{ "Book" output.bibitem }

FUNCTION {booklet}{ "Booklet" output.bibitem }

FUNCTION {inbook}{ "InBook" output.bibitem }

FUNCTION {incollection}{ "InCollection" output.bibitem }

FUNCTION {inproceedings}{ "InProceedings" output.bibitem }

FUNCTION {conference}{ "Conference" output.bibitem }

FUNCTION {manual}{ "Manual" output.bibitem }

FUNCTION {mastersthesis}{ "MastersThesis" output.bibitem }

FUNCTION {misc}{ "Misc" output.bibitem }

FUNCTION {phdthesis}{ "PhdThesis" output.bibitem }

FUNCTION {proceedings}{ "Proceedings" output.bibitem }

FUNCTION {techreport}{ "TechReport" output.bibitem }

FUNCTION {unpublished}{ "Unpublished" output.bibitem }

FUNCTION {default.type} { misc }

MACRO {jan} {"January"}

MACRO {feb} {"February"}

MACRO {mar} {"March"}

MACRO {apr} {"April"}

MACRO {may} {"May"}

MACRO {jun} {"June"}

MACRO {jul} {"July"}

MACRO {aug} {"August"}

MACRO {sep} {"September"}

MACRO {oct} {"October"}

MACRO {nov} {"November"}

MACRO {dec} {"December"}

MACRO {acmcs} {"ACM Computing Surveys"}

MACRO {acta} {"Acta Informatica"}

MACRO {cacm} {"Communications of the ACM"}

MACRO {ibmjrd} {"IBM Journal of Research and Development"}

MACRO {ibmsj} {"IBM Systems Journal"}

MACRO {ieeese} {"IEEE Transactions on Software Engineering"}

MACRO {ieeetc} {"IEEE Transactions on Computers"}

MACRO {ieeetcad}
 {"IEEE Transactions on Computer-Aided Design of Integrated Circuits"}

MACRO {ipl} {"Information Processing Letters"}

MACRO {jacm} {"Journal of the ACM"}

MACRO {jcss} {"Journal of Computer and System Sciences"}

MACRO {scp} {"Science of Computer Programming"}

MACRO {sicomp} {"SIAM Journal on Computing"}

MACRO {tocs} {"ACM Transactions on Computer Systems"}

MACRO {tods} {"ACM Transactions on Database Systems"}

MACRO {tog} {"ACM Transactions on Graphics"}

MACRO {toms} {"ACM Transactions on Mathematical Software"}

MACRO {toois} {"ACM Transactions on Office Information Systems"}

MACRO {toplas} {"ACM Transactions on Programming Languages and Systems"}

MACRO {tcs} {"Theoretical Computer Science"}

READ

% DFK use cite$ for sort key
FUNCTION {presort}
{
  cite$
  #1 entry.max$ substring$
  'sort.key$ :=
}

ITERATE {presort}

SORT

FUNCTION {begin.bib}
{
  "% BibTeX bibliography file" write$ 
  newline$ 
  newline$
}

EXECUTE {begin.bib}

ITERATE {call.type$}