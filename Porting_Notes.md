# 3photons-julia experiment report

## Introduction

### Context

The world of numerical computing is a lively place. All year long, new
algorithms are devised, new libraries are released, and new programming
languages and paradigms promise to change the way we think and code for the
better, be it towards more efficiency, more portability, or less mistakes.

Diving through this intense intellectual activity, however, can be daunting for
the practicioner, and especially for those of us for whom the craft of computing
is not a primary activity but more of a secondary tool, a mean towards an end
that should never become an end in and of itself.

Therefore, it is part of the duty of technical personnels who are most willing
to expend time and energy into the study of these evolutions to report on their
findings, and to provide concrete advice and recommendations to their less
technically-minded colleagues regarding which of these new approaches seem most
promising in the long term and worthy of further investment from their side.

Which is how I ended up evaluating the viability of Rust for numerical computing
previously, and am now evaluating that of Julia, through the interesting
exercise of porting a simple but non-trivial Monte Carlo simulation to this new
programming language and taking notes in the process.

### Guinea pig

"3photons" is a Monte Carlo particle physics simulation. In a nutshell, it
generates many random occurences of a certain physical process of interest, then
uses these random data points to approximately integrate some physical
quantities which characterize that process, such as its cross-section
(like a probability of occurence, but without experimental setup effects).

This simulation, which was initially developed by my colleague Vincent Lafage as
part of his PhD work, has several characteristics that make it interesting as a
technology testbed for numerical computing tools :

- At a couple thousands of lines of code, it is simple enough to be ported to a
  new technology in a couple of weeks while doing other things concurrently, yet
  complex enough to start stressing some advanced aspects of programming
  languages such as support for keeping larger codebases tidy or cleanly
  expressing abstract mathematical concepts in code.
- It stresses several common aspects of numerical codes, including random
  number generation, complex number support, small-matrix linear algebra,
  readability of complex mathematical expressions, and controlled text I/O.
- It has been ported to many different programming languages (Fortran 77 and 90,
  Ada, C, C++, Go, Rust, and now Julia) and paradigms (OpenMP, OpenCL, HPX)
  before, so there is ample experience about what the typical areas of
  differentiation to watch out for are.
- Much energy has been expended into optimizing the program over the course of
  these ports, so the core algorithm is known to exhibit very good CPU cache
  locality and parallel scaling properties. Any failure to uphold this
  expectation can therefore be safely attributed to either the quality of the
  port or a problem of the technology stack being evaluated.

My goal, through this 3photons porting exercise, will be to evaluate how Julia
seems to behave on a non-toy benchmark, in areas including...

- Run-time performance and "performance ergonomics", i.e. how easy it is to
  write code that runs fast, and whether performance is easy to understand and
  optimize or to the contrary code contortions are needed for good performance
- Small-scale code clarity, e.g. degree of math expression clutter induced by
  syntax noise like explicit int->float conversions, availability of
  prefix-notation operators which are clearer to scientists...
- Large-scale code organization abilities, e.g. module system ergonomics and
  ability to express domain concepts as code abstractions...
- Ergonomics of writing code that is both correct and hard to use incorrectly
- Built-in and external lib support for common numerical computing scenarios
- Quality of documentation and learning materials, overall ease of learning

My starting point and main point of comparison will be the Rust version of
3photons at https://github.com/HadrienG2/3photons-rust because...

- I wrote this one and therefore have the deepest understanding of it
- It is one of the versions that has received the most validation and
  performance optimization care, which sets a good bar for the Julia version
- It is a good occasion to compare the current state of the art in
  performance-minded "static" and "dynamic" languages.

### Personal bias disclaimers

As a human being, I make mistakes. I love to work with tools that reliably catch
my mistakes before they go into production, and I don't mind spending a bit of
time to learn how to use these extra tools in exchange. In a programming
context, this means that I am very much a fan of advanced static type systems,
and more dynamic systems which either catch my mistakes at runtime or don't
catch them at all will start with a negative apriori from my side.

As a software engineer, I am ready to dedicate significant time to learning a
new technology, so I am less sensitive to "how quickly can I learn X" concerns
than many others will.

I am also very much used to file-based programming workflows (together with
proper code/data separation and version control), and averse to interactive
workflows like REPLs or notebooks. So I will only use the latter for quick
checks and where necessity dictates (as when e.g. installing packages in Julia),
which I think does not do full justice to the great work that the Julia
community has done in this area. My workflow may instead stress unusual and
underdeveloped aspects of the language and its implementation.

As a performance engineer, I tend to have a very severe opinion of any
programming language which doesn't perform as well as has been proven possible
in another language, with double negative points if the implementation gets in
my way as I try to understand the origin of the difference and bridge the gap
through optimization. This thinking is part of what got me interested in Julia
over the current Python statu quo, but it can also work in Julia's disfavor when
some of its dynamic features prevent it from reaching the same level of
automatic optimization and performance ergonomics that its statically typed and
AoT-compiled cousins enjoy.

To speed up porting and ease comparisons by colleagues, I closely mimicked the
Rust code's organization on the Julia side, which may at times have come at the
expense of following idiomatic Julia style and constructs. Any opinion and
advice on this front is much welcome, although I will be a little biased against
following nontrivial code organization suggestions for the aforementioned sake
of keeping the two versions as comparable as possible.

Finally, although I personally think that small matrix algebra is an important
feature of a numerical computing language, I must admit that this has proven
time and time again to be a niche feature which is only supported via relative
obscure libraries, such as `StaticArrays.jl` in Julia's case. As a result, I may
occasionally misattribute Julia blame for something which is squarely the fault
of StaticArrays, even though I took some care to put credit and blame where it's
due on this front.

Now that his is out of the way, without more ado, let's go into my notes.


## Documentation

I found the documentation at https://docs.julialang.org/ reasonably complete but
very hard to use, both for learning purposes and as a reference, because :

- It frequently switches between a verbose tutorial style and a concise
  reference style, in spite of these two forms of documentation having very
  different requirements and use cases.
- It lacks structure. It is not unusual for a single section to be several pages
  long without any subsection. This makes it hard to interrupt reading and
  resume it later, or to find back an information that you read in the past.
    * Reference pages are especially bad in this respect, there are often just
      very long lists of methods with links to the source code, with low visual
      separation between individual method entries.
- It lacks consistency. Each page feels like it was written by its own author in
  its own special way, without an overarching editorial policy to make page
  structure predictable and ease information lookup.
    * In addition, I find the very uneven level of detail depending on the
      subject matter very problematic. It could be better to have, for example,
      a beginner tutorial where everything is lightly touched upon, followed by
      an advanced tutorial where topics that need it are explored in detail.
- The little structure that exists is neither announced at the beginning of a
  page, nor easily navigable through hyperlinks in the sidebar, which only go
  through the first level of headings (which can span many pages of content).

The search bar helps one work around some of these structural deficiencies, but
only up to a point:

- Every search query takes uncomfortably long by modern search engine standards.
- As with all search-based queries, if you don't remember a keyword that is
  clearly associated with the concept you have in mind, you are out of luck.
- Some search terms can yield many different results, and the search result page
  does not do a very good job at explaining why a specific query matches in a
  specific area. Some small quotes with search term highlighting, as e.g.
  Google and DuckDuckGo do, could help greatly here.

Third-party documentation sources like the ["Introducing Julia" wikibook](
https://en.wikibooks.org/wiki/Introducing_Julia) were usually not much better
organized than the official manual, and occasionally out of date, so I didn't
expend too much time into perusing them unless I was looking for a specific
topic that wasn't clearly covered by the Julia manual.


## Small-scale syntax

### Expression clarity

The Julia language shines where it is most strongly expected to deliver, namely
at making math expressions feel as much like pseudocode as possible, without
making the language grammar too ambiguous in the process.

Here are some language and community decisions that contribute to this feeling:

- The local culture of using Unicode identifiers like greek letters and math
  symbols goes a long way towards making Julia code look like math. It does do
  so, however, at the cost of making code harder to edit (you need either
  copy-paste, a charmap, or a specialized code editor/REPL), so I would only
  recommend this as a late code post-processing step before publication.
- The abstract hierarchy of math types makes it very convenient to write code
  that is a bit generic, but not too much, e.g. when narrowing down when
  dispatch on a method should fire or not.
- The conversions and promotion rules are remarkably well thought-out, and most
  of the time typing the equivalent of a math expression will Just Work, without
  any need for explicit conversions, markers of float-ness (e.g. typing `3. + x`
  instead of `3 + x`) or C-style broken promotion rules that fire too late after
  some information has been lost (e.g. `1 / 3 * 3` returning `0`).
- The numerical prefix syntax `3x` initially felt like a gimmick that makes the
  language grammar needlessly ambiguous (e.g. can you remember in a fraction of
  a second what `7/8x` does? And if the goal is to look like math, why do I need
  to write `7x*y` over `7xy`?). But after using it, I must admit that it does
  help reduce syntax noise and expression width in many situations.

All these design touches contribute to making numerical expressions very concise
and light on syntax noise in the common case, which makes the mathematical
intent behind them clearer. Again, great work here.

### Literal typing

After a few years of using civilized programming languages that treat number
literals as members of the abstract mathematical set of integers and reals, as
opposed to as machine numbers of a specific signedness and width, it pains me to
once again use a programming language where passing the literal `1` to a piece
of code that expects an `UInt` or passing the literal `4.2` to a piece of code
that expects a `Float32` leads to a type error.

However, I am well aware that this problem is most likely caused by Julia's
decision to base everything on multiple dispatch. Like all forms of function
overloading, it makes input parameter type inference undecidable. So, language
design tradeoffs as usual, I guess...

### Array indexing

I am deeply saddened by Julia's decision to favor 1-based array indexing and
right-inclusive ranges by default.

Divergences across programming languages on array indexing and range conventions
have caused innumerable off-by-one problems as programmers had to move from one
language to another, or worse, work on multi-language projects. For a long
while, it felt like the programming community was finally in the process of
standardizing on zero-based indexing and right-exclusive ranges. It pains me
greatly to see a new programming language throw more oil on this old fire that
was almost extinct, and thus keep this mess alive for future generations.

And no, allowing user-defined indexing conventions is not the answer. The only
thing this language feature achieves is to make array indexing code almost
impossible to write correctly, as devs will work with 1-based arrays 99.9% of
the time and never remember to test their code against arrays that use other
indexing bases.

Still on the matter of indexing, it also feels strange that arrays are indexed
with `Int`, and not by `UInt`, given that Julia does not use negative indexing
like e.g. Python does (which is for the best, don't get me wrong).

It seems to me that this is a sad byproduct of another Julia design decision,
which I already lamented above, that made integer literals typed instead of
inferring their types based on usage like any modern Hindley-Milner type
inference engine is able to.

### Docstrings

The Julia docstring syntax holds the dubious distinction of feeling even more
like a bad parser hack to me than the Python one. But that is arguably personal
taste. What does not boil down to taste, however, is that Julia's docstring
syntax inherits the same usability problems that Python docstrings have been
facing forever:

- Unless your code editor's Julia parser is sophisticated enough to distinguish
  docstrings from other strings, which is far from a trivial task, their visual
  highlighting will be close to that of string expressions in method
  implementations, and far away from that of comments. This gives the wrong
  ergonomics message as docstrings are code documentation, not data.
- Docstrings tend to be subjected to all kinds of special internal conventions
  that do not concern regular strings. It is very hard for code editors to help
  you visualize and follow these extra semantics, since they must resolve the
  aforementioned string/docstring disambiguation problem before doing that.

### Whitespace-delimited grammar considered harmful

Over time, I've gradually come to dislike use of whitespace over more explicit
delimiters (commas, semicolons, keywords...) in programming language grammars.

At first, it sounded like a net ergonomic improvement brought by modern parsing
techniques. Finally, beginners would be freed from the insufferable burden of
remembering to type a semicolon at the end of every line of code!

But then everyday experience with using languages that made this design choice
taught me that whitespace as a delimiter brings with it many papercuts that
cause everyday pain to practicioners, a pain which feels greater than the
benefit of eliminating the one-time cost of learning to correctly spell
expression delimiters as a beginner. Let me spell out some of its problems.

---

First, whitespace as a function/macro argument separator makes any nontrivial
expression ambiguous. For example, in `@something a b + c`, it's not immediately
obvious to a reader with a normal mind whether `@something` takes two arguments
`(a, (b + c))` or four arguments `(a, b, +, c)`.

I know that the grammar is not ambiguous in a literal sense. In the end,
operator precedence rules will decide what the implementation will do. What I'm
saying is that whitespace as an argument separator makes it trivial to write
code that is hard for human to read, unlike more explicit separators like
parentheses and commas.

---

Second, whitespace as a block delimiter tends to facilitate code typos. For
example, in Python or Nim...

```python
if order_given:
    load_missile()
launch_missile()
```

...may not do what the author expected, and may be hard to spot in the
messy identation style of inexperienced beginners (for the benefit of whom
we're enduring all this extra parser complexity to begin with!).

I believe that Julia is immune to this problem, however, by virtue of
strongly favoring line feeds after conditionals in control flow expressions and
mandating explicit `end` block terminators.

---

Finally, whitespace as a statement delimiter forces complex expressions to be
broken on awkward boundaries. For example, consider the following expression,
assuming that `a`, `b`, `c` etc. are actually big math expressions for which
a line break is warranted...

```julia
(a +
 b -
 c +
 d) /
(e +
 f -
 g +
 h)
```

...and now, compare with this style:

```c
(a
 + b
 - c
 + d)
/ (e
   + f
   - g
   + h)
```
The latter style is, in my biased opinion, much easier to read than the former.
But unfortunately, enforcing whitespace as a statement delimiter makes this
style illegal under the language's grammar...

### Macro syntax inconsistency

It is unfortunate that the syntaxes for declaring and using a macro are so
different from each other (`macro something(foo, bar, baz)` vs
`@macro foo bar baz`. I think there's a small lost opportunity for syntax
consistency here...


## Large-scale code organization

### Code organization tools

I am not very impressed by the mechanisms provided by Julia to organize large
codebases which cannot reasonably fit in a single file, nor by the way in which
these constructs are introduced in the official documentation :

- While it is occasionally useful to have a textual inclusion mechanism that's
  decoupled from the module system, it is a bit of a specialized tool. For most
  purposes, 1 module = 1 file is good enough, and more ergonomic as...
    * It reduces redundant syntax (no need for 1 piece of textual inclusion
      syntax somewhere *and* 1 piece of module definition syntax somewhere else)
    * When combined with explicit imports, it makes it easy to know in which
      file a certain construct is defined.
- The current setup makes it hard to track inter-module dependencies. We cannot
  just `include()` a module in every other module that uses it, because that
  will lead to duplicate declarations, which will in turn lead the Julia
  implementation to rightfully complain.
    * In many respect, this situation is very similar to the C header situation.
      Hacks similar to the `#ifndef BLAH/#define BLAH/#endif` triad can probably
      be applied in Julia too. But one would have thought that 50 years after C
      was designed, programming language designers should have finally learned
      that textual inclusion makes a very poor dependency tracking system.
    * Further, the C header system has several other problems which Julia is
      prone to experience, such as "implicit" inclusion of a module by a
      virtue of this module being ordered earlier in the final `#include` list.
- The package system documentation, which is mostly featured in the somewhat
  unfortunately named [Code Loading](
  https://docs.julialang.org/en/v1/manual/code-loading/#Federation-of-packages-1)
  section of the Julia manual, is diluted by heaps of complex decentralized
  system discussion, which makes it feel like a very advanced topic. So I am not
  convinced by the argument, made in e.g. [this blog post](
  https://white.ucc.asn.au/2020/02/09/whycompositionaljulia.html), that
  packages are the first code organization tool that Julia devs should go for.
    * It is not enough to tell people to go for abstractions like
      [PkgTemplates](https://github.com/invenia/PkgTemplates.jl/) and ignore the
      detail of what's going on either, because if they don't have a robust
      understanding of the underlying system, they will not be able to avert
      default behavior from the abstraction that don't do what they want, or to
      understand problems with the base layers when they do occur.
    * Further, the argument made by this author that every programming language
      design problem can be resolved through developer communication and
      sufficient mutually-agreed-upon conventions has actually been made many
      times before, in the context of many other then-young programming
      languages, and systematically shown by experience not to scale well to
      larger programmer communities.
- The `import` vs `using` distinction is very subtle and easy to misunderstand,
  yet getting a complete picture of it requires jumping through various pages of
  manual and reference documentation, plus reading through some external
  documentation like the Julia wikibook.
    * Given how crucial a difference it makes to multiple dispatch, which is a
      central feature of Julia, I feel like this should be explained better, not
      just left to an obscure FAQ entry that most beginners will miss.
- The `export` keyword is introduced misleadingly in the official documentation.
  Said documentation strongly gives the wrong impression that modules can act as
  a privacy boundary, whereas further investigation reveals that `export` only
  controls what is made visible by blanket `using` statements, and not what can
  be accessed via more explicit `Module.stuff` syntax.

### `export` and privacy barriers

As I was deeply disappointed by Julia's `export` mechanism, let me complain a
little more about it.

I consider Julia's lack of serious privacy barriers problematic because it
hampers package evolution. In programming languages with a proper privacy system
in place, you can confidently do things like add new private members to structs
(e.g. for caching purposes), as you know it won't break your API's clients.

In Julia, you cannot do this, because there isn't even a way to prevent users of
your package from calling a struct's default constructor, which will break if
you add one new struct member.

Even if we consider `export` as a gentleman's agreement ("here's what you
_should_ be using, I can't guarantee backcompat if you touch other things"), it
is still less convenient in this role than `public:` sections in C++-style class
declarations or `pub` qualifiers in Rust, because it forces you to either have
lots of `export` statements or write the information about what's being exported
by a module far away from the point where the thing that's being exported is
being declared.

### Sensitivity to declaration order

Still on the topic of Julia feeling like a regression compared to other
contemporary programming languages in some places, I should also mention the
language's awkward sensitivity to entity declaration order, which feels like a
knockback into a long-forgotten punchcard era where compiler authors couldn't
afford processing programs in more than one pass.

- Once one has gotten used to it, the ability to write code modules by putting
  important user-facing APIs first, and less relevant implementation details
  second, gives code a degree of clarity that is hard to give up on.
- Some kind of programs, such as those based on mutually recursive functions
  (think a graph traversal engine with one function for traversing nodes and one
  function for traversing edges), are unpleasantly hard to write in programming
  languages that are sensitive to function declaration order.
    * ...and unlike in C/++, Julia doesn't seem to provide a way to declare a
      function without defining it, which is a common way to work around this
      issue at the cost of some code redundancy.

I do not know if the decision of having module processing be sensitive to
declaration order was motivated by ergonomics consideration (keeping module
processing consistent with direct REPL input) or implementation convenience, but
I wish this decision could be revisited someday.

### Multiple dispatch and method namespacing

A discussion of Julia's strong and weak design points would not be exhaustive
without a discussion of multiple dispatch, namespacing, and how these two
programming language features interact with another. Here's my take:

- `3photons-julia` does not do anything fancy enough to benefit from internal
  use of multiple dispatch (except for a tiny occurence of specialized text
  output), so I can only comment on how I saw other people use it.
- Multiple dispatch-based APIs tend to use a mixture of overloading and duck
  typing. Experience with similar combinations in C++ shows that this
  combination, while powerful, can be hard for programmers to reliably reason
  about. This is especially true when the method overload set is open for
  extension and can therefore change depending on which files were `include()`d
  beforehand, as is the case in both C++ and Julia.
- While much virtual ink has been spilled about how multiple dispatch and module
  namespaces interact (as the former tends to favor a "flat" API style that
  makes the latter feel out of space, resulting in frequent occurence of blanket
  `using Something` in Julia), I felt that what affected me most in my personal
  code is the loss of _class_ namespaces, i.e. methods being attached to types:
    * It means that when using explicit imports, you must import many things
      instead of ones. This makes such imports harder to use, more verbose.
    * It means that method names must often be made more verbose, because you
      cannot count on the receiver to act as a context disambiguator as in
      single-dispatch language. For example, `event_cut.keep(event)` is easily
      understandable to any HEP Monte Carlo expert, while
      `keep(event_cut, event)` feels insufficiently clear and one feels the need
      for a more verbose `keep_event(event_cut, event)`.
    * It means that methods that implement interactions between two types of
      object have _even less_ of a natural place to go to than in
      single-dispatch languages. Should they be defined in the module where the
      first object is defined? In the module where the second object is defined?
      In a third module, even though that will make imports less obvious?

### Macro hygiene

Here is a `StaticArrays` user issue, which feels like it _might_ actually
originate from a Julia design issue, but I may be mistaken: why is it than when
I write `using StaticArrays: @SVector` macro, and then try to use the macro, I
get an error that tells me to bring the `SVector` struct into scope as well?
This feels like an undesirable macro implementation detail leakage to me.

The issue might originate from the design of Julia, rather than from the
implementation of the `StaticArrays` package, because I don't think I've seen
Julia support for generating import statements in functions, which I believe
the `@SVector` macro would need to do in order to save me from the trouble of
adding an extra import myself.

From what I understand, it's a combination of several Julia design decisions...

- Needing to type in an explicit using/import statement to use module elements.
- No support for using/import statements in functions.
- Limited macro engine support for macro declaration site hygiene.

...but a comment from a StaticArrays dev or Julia designer that confirms or
contradicts my understanding here would be much appreciated.

### Module vs content naming

It is somewhat sad that in Julia, modules and types are expected to follow the
same CamelCase naming convention, but a module cannot bear the same name as an
inner struct (as it will lead a constant redefinition error).

When a module is 100% dedicated to implementing a single complex type, as is
often the case, it would be tempting to name it in the same way as the type that
it is implementing, instead of coming up with an artificially different name
that users will also need to remember.

In `3photons-julia`, examples of that include `EvData.Event` and
`EvCut.EventCut`, where it would have been more tempting to name the former
module `Event` and the latter module `EventCut`.


## Type system and correctness checks

### Named struct constructor syntax

When a Julia struct is defined, a constructor function that takes all fields in
declaration is automatically generated. It would be nice if a variant of this
constructor function that uses named argument syntax, with the struct's field
names as named arguments, could also be automatically generated :

- Usage of named arguments clarifies which constructor argument maps to which
  struct fields. This is helpful when a struct's definition is quite remote from
  usage of its constructor, or when a struct has many members.
- Usage of named arguments makes the codebase more robust in the face of future
  evolution, where seemingly innocent refactorings that put struct fields in
  a more logical order (or in alphabetical order) can break callers.
- Absence of named arguments greatly increase the odds of a common code typo
  where struct fields are listed in the wrong order during constructor call.

To summarize, struct constructors with named fields have some advantages, which
for some people are worth the trouble of typing field names in constructor
calls. It would be nice if this were natively supported by the language.

### Functional-first design tradeoffs

I appreciate the care that has been taken to make the functional programming
style ergonomic and performant in Julia, as I feel it is often the right default
for numerics code. However, this should not come at the expense of making
imperative style harder to write and less performant, as currently happens due
to a combination of...

- Structs being immutable by default (as opposed to being mutable but normally
  passed by read-only _references_ or _bindings_, as in e.g. Rust).
- Mutable structs being allocated on the heap unless the compiler manages to
  avoid it, which has serious performance implications.
- Structs being often defined in two closely related mutable and immutable
  versions, like e.g. `MVector` and `SVector` in StaticArrays, which leads to
  much code/API duplication and, again, longer explicit imports.
- Iterators being quite unpleasant to use in an imperative way, e.g. to mutate
  an array that is being iterated upon. Even explicitly calling the iterator
  protocol without for loop sugar is quite an unpleasant experience in Julia
  compared to most other programming languages with built-in iterator support.

### Type system advantages of dynamism

For all I could complain about how Julia's dynamic nature makes performance
optimization painful, I must admit that it does have ergonomic advantages:

- Since the compiler is always around, instantiating generic code based on
  run-time configuration is not a problem. No need for pre-instantiations and
  the like, you just call the thing and it magically works.
- Not only can types be manipulated like values as in languages with
  introspection, the reverse is also true. Generic code (parametric code in
  Julia jargon) accepts values as parameters just fine, with much less
  complications than in many other languages that accept values as generics
  parameters. For the most part, you put a `Val` on it and it just works. This
  makes things like fixed-sized arrays much more convenient to implement & use.

### Inefficient static error reporting

When I write Julia code, I miss the ability of statically typed programming
languages to _exhaustively_ check for a whole category of error in a single
compiler run. Debugging freshly written Julia code without a REPL is a
frustrating experience by comparison:

- Start the program.
- Wait 10s for the JIT to be done messing around.
- Get a single `UndefVarError` about a typo somewhere.
- Fix the typo, restart the program.
- Wait 10s for the JIT to be done messing around.
- Get another `UndefVarError` about another typo elsewhere.
- Fix the typo, restart the program.
- Wait 10s for the JIT to be done messing around...

...I think you get the point. I don't mind Julia supporting people who like
REPLs even if I personally don't, but I take issue with how badly Julia punishes
those who are not using the REPL with those huge turnaround latency and
error-by-error issue reporting.

### Lack of exhaustive switches

As someone whose main programming languages are C++ and Rust, I miss switch
statements in general, and ML-style exhaustive pattern matching in particular.

It feels wrong to type in long elif sequences, and have them terminated with an
`else` that throws an exception because without that I have no way to know when
I messed up consistency between the enum declaration and my switch statement.

### Callable typing constraints?

I wonder if one can express precise typing constraints on callables in Julia,
like e.g. `impl Fn(usize, &str) -> bool` does in Rust. So far, all I found was
the `Core.Function` abstract type, which as far as I can see doesn't allow
specifying what the expected arguments and results of the function are.

In addition to hampering use of type checking as a program consistency
assertion (which I'm well aware is not a very high-priority use case of type
annotations in Julia), this could also prevent some advanced uses of multiple
dispatch where different versions of a higher-order function are called
depending on what kind of function is passed as a parameter.

### Type instability inviting more static typing

I would like to point out one possible future benefit of annotating function
input and return types, besides making static typing fans happy and enabling
multiple dispatch: if the compiler supported this use case well, it would enable
catching many type instability performance issues, without going through the
hassle of manually running `@code_warntype` through every function of the hot
part of your codebase or trying in vain to make sense of `Traceur`'s output.

Due to the peculiarities of mathematics in general and floating point in
particular, it's not unusal in numerics to have code that looks like this:

```julia
result = if is_normal_case(x, y)
    do_work(x, y)
else
    SPECIAL_VALUE
end
```

That pattern happens for many reasons:

- Fast-path optimizations, in special cases where a simpler expression of the
  desired mathematical result is known.
- Regularization of numerically unstable operations, such as division of a small
  number by another small number.
- Genuine discontinuity in the mathematical behavior of the computation under
  study, as in e.g. step functions or absolute values.

But in Julia it's all too easy to get it wrong, as it's a situation where the
language's extremely well thought-out integer-to-float promotion rules, which
make use of integers in place of floats a non-issue most of the time, come back
to bite us. For example, this is type-unstable and slow:

```julia
function ramp(x::Float32)
    if x > 0
        x  # Returns a Float32
    else
        0  # Returns an Int
    end
end
```

The compiler, however, should be able to report the problem automagically if we
were more explicit about the kind of output that we expect from the `ramp`
function:

```julia
function ramp(x::Float32)::Float32
    if x > 0
        x  # Returns a Float32, OK
    else
        0  # ERROR: Returns an Int, but a Float32 was expected
    end
end
```

However, this currently produces neither a compile-time nor a run-time error.
Does anyone agree with me that it probably should?

### Need for error reporting work?

This above was an example of a pattern which I observed more frequently in
Julia, where obviously erronerous code would not only parse without a warning,
but execute without a runtime error. Another example was `export`ing nonexistent
entities. Unfortunately, I didn't take rigorous notes of those, as I was
generally too quick to fix the code in anger and move on to other things, so I
cannot provide a more exhaustive list.

But I got the impression that Julia's error detection and reporting story is a
bit spotty right now and this topic could use more rigor in the future.

### Thoughts on conditional compilation

I'm a bit on the fence regarding Julia's support for conditional compilation, or
more precisely its lack thereof.

On one hand, it makes sense in the context of a dynamic language where the
boundary between run-time and compile-time is fuzzy by design anyway. Which, as
I mentioned earlier, does make some things easier than "static" designs where
there is a strong distinction between compile-time and run-time concerns.

On the other hand, it means that "run-time" function calls must be polluted with
"configuration" concerns, such as the choice of floating-point type in use or
whether a given simulation feature (like multi-threading) is enabled.
Compile-time configuration are one of those situations where I personally think
that global objects are the least of available evils.

In principle, this global object design could be replicated by implementing
"compile-time" settings as a separate file in the source full of `const`s that
users are expected to edit. But that feels decidedly more troublesome for users
trying to experiment with the various compile-time settings than being able to
simply add some flags to the julia command-line invocation.


## Libraries and APIs

### Math coverage from the standard Julia library

As expected from a modern numerics language, Julia does a good job of providing
every math operation I could possibly need in this project, even though I needed
to use `StaticArrays` for performance reasons.

The ability to use every math operator in prefix form without extra tricks is
a strong point of Julia's unusual choice not to use OO-style `receiver.method()`
syntax, even though that choice was most likely rather motivated by dynamic
dispatch considerations.

However, this lack of OO module structuration does tend to negatively impact
documentation, since many modules feel like a huge unsorted heap of method. More
tools should probably be provided to package authors so that they can make their
automatically generated docs easier to navigate by structuring their method list
better.

### Text output formatting control

Custom text output formatting was a surprisingly low point of this project.

- Apparently, `@printf` is the only available tool for this purpose in the
  standard library. But its C-ish formatting DSL feels totally out of place in
  Julia code, and it does not even provide full feature compatibility with
  advanced features of the C printf function, such as dynamic width and number
  precision specifiers.
    * This particular absence would make it a headache to write a program that
      can switch between various floating-point precisions and prints the right
      number of significant digits at the end, which is one future extension
      that `3photons-julia` would need for feature-parity with the Rust version.
- I saw the [Formatting](https://github.com/JuliaIO/Formatting.jl) package
  being sometimes recommended as an alternative, but unfortunately it did not
  impress. Its documentation can be summarized as "just type in a valid Python
  format string and it should work", and that's clearly false:
    * `g` formatting is not supported. That's documented, but a big issue in
      scientific computations.
    * Padding is not supported. That's not documented, and a big issue when
      printing tabular output.
    * Overall, I'm not sold at all on the core idea of cloning Python's DSL. Now
      I need to learn both said DSL and the ways `Formatting.jl` deviates from
      it, and the compatibility gap is likely to increase as both Python and
      `Formatting.jl` continue to evolve.

### `println()`'s leading-stream optional input

I'm not sure if I'm sold on println's design decision of handling the decision
of printing to `stdout` or another IO stream via an optional _head_ argument:

- It's at odds with the usual convention of putting optional arguments at the
  tail of the argument list.
- Using the same function name for `stdout` and file output makes it somewhat
  too easy to accidentally print to `stdout` by forgetting the stream input.
    * On the other hand, these sorts of errors should be relatively easy to
      catch during testing most of the time, so maybe it isn't too bad.

### StaticArrays thoughts

I was intrigued to see how Julia would handle fixed-sized matrices, as
experience with other programming languages (Eigen in C++, nalgebra in Rust...)
told me that implementing fixed-sized matrices in a language that doesn't have
built-in support for it can be trickier than it sounds.

StaticArrays is what I identified to be the current dominant Julia library for
this problem, and overall I was mostly satisfied with it. It does a good job of
"feeling native" overall, while delivering on the expected goal of avoiding
heap allocations in the MC simulation kernels (at a fraction of a µs per
iteration, the allocation loop cannot affort much of it).

As in other areas of Julia, however, I was a bit saddened by the amount of
"performance traps" that StaticArrays exposed:

- It's a shame that slicing on a range of indices that's known at compile-time
  (e.g. `arr[2:4, 3:4]`) cannot be optimized out as well as
  `arr[coord, :]`-style slicing, and that weird constructs like
  `arr[SVector{3}(2:4), SVector{2}(3:4)]` are apparently needed to get good
  performance out of such slices.
- It's sad that the compiler doesn't manager to optimize out the heap
  allocations of mutable array types like `MVector` except in trivial cases
  (e.g. a row swap is all it takes to break optimizations, as does an hcat/vcat
  of `MMatrix` although the result is immediately converted to `SMatrix`...)
- It's ugly that number of matrix elements must be specified _in addition to_
  number of matrix rows and matrix columns in order to avoid type instability
  when StaticArrays are used as struct members.
- The docs could provide more advice on these performance traps and the ways in
  which one can avoid them.

### The mutation_bang_suffix! API convention

It is also not clear to me if the bang!-suffixed coding convention for methods
that mutate one of their arguments is worth the trouble, especially as [the
recommended method argument ordering](
https://docs.julialang.org/en/v1/manual/style-guide/#Write-functions-with-argument-ordering-similar-to-Julia-Base-1)...

- Does not guarantee that the input(s) being mutated will always be in a
  well-identifiable position, such as "always first".
- Does not clarify if multiple inputs are being mutated, and if so how many.

On this front, I guess I prefer the way idiomatic Rust APIs effectively force
the caller to explicitly annotate function arguments with are being mutated with
`&mut`. If you are going to use loud syntax for mutation (which, in and of
itself, not everyone will be sold on), a Rust-style mutated _argument_ marker
seems like a better way to do it.


## Performance ergonomics

### Performance vs code organization

On the topic of language design decision unpleasantly affecting the daily
development experience, I would like to give a mention to the many ways in which
code structure interacts with run-time performance in Julia, as abundantly
discussed in the [Performance Tips](
https://docs.julialang.org/en/v1/manual/performance-tips) section of the Julia
manual:

- Global variables basically break static type inference, and `const` does not
  always recover from all the damage.
- Closure captures basically face the same problem plus extra boxing overhead,
  which makes closures hard to use in cases where the compiler's optimizer
  doesn't manage to eliminate this overhead.
- The fact that struct mutability affects the compiler's stack/heap data
  placement decisions feels like an undue burden on code authors, even when one
  knows the dynamic typing rationale for it. Not all code is fit for being
  written in a functional style.
- The effect of type annotations on the compilation process varies depending on
  context. On structs fully specific type annotations are almost mandatory for
  performance, while on function declarations they are only useful as assertions
  and markers of multiple dispatch. This yields idiomatic Julia code to have an
  inconsistent coding style where some types are annotated and others aren't.
- Function boundaries are a bit too special in the eye of the Julia compiler,
  and I could easily see cases where developers are forced to distort their code
  into unnatural shape in order to achieve optimal run-time performance, much
  more than is the norm in AoT-compiled languages:
    * You may need to extract code into a separate method in order to get more
      specialized (and therefore more performance) code.
    * You may need to manually inline a method into another because the
      compiler's current inlining heuristics are very conservative and the
      `@inline` annotation has strange side-effects like causing allocations.
    * As mentioned, the performance ergonomics of closures are rather low.
- The fact that modules act as a precompilation boundary leads to [scary and
  unexpected cognitive burden](
  https://docs.julialang.org/en/v1/manual/modules/#Module-initialization-and-precompilation-1)
  as soon as one has to go beyond the common case where the difference between
  precompiled and non-precompiled code is not observable.

In this respect, I feel like Julia fails to live up on its promise to bring the
performance ergonomics of dynamic typing close to those of a statically typed
language through its use of JIT-compilation techniques.

Not only does one often need to write code that "feels" statically typed to get
good performance, the mere _possibility_ of using dynamic typing restricts
compiler code optimizations and leads to counterintuitive implementation
behavior even when that possibility is not being used.

Maybe providing more static typing annotations (e.g. attesting that a global or
captured variable will always have the same type at its point of declaration
instead of at its point of use) could help a bit here?

### JIT performance

I must join the choir of complaints lamenting the state of Julia JITting times.

Even on this small-ish 2kLOC project without much generic code, I face total
running times of 15.2s of which the actual (pre-warmed) computational kernel
only accounts for 5.8s, and this 1:2 "signal-to-JIT ratio" feels a bit much.

Add to that that any performance analysis tool (even `@time`!) contributes quite
a bit to running times in a fashion that measurably biases the observed
performance profile, and hopefully you can see why I am unhappy.

I wish the Julia team would replicate what the Rust team did when facing build
time issues, and start keeping around caches of JITted codes across program
runs. I have the impression that the `--compiled-modules=yes` flag is intended
to be a step in that direction, but it does not seem to have much of an effect
on my code base as of Julia 1.3.1.

Maybe the cache is too coarse-grained and too readily invalidated by minor
program or environment changes? Or maybe I'm just misunderstanding what this
module incremental precompilation business is supposed to be doing, and this
cache is only used in interactive scenarios like REPL & notebook use?

From an outsider perspective, this lack of "offline" JITted code cache just
feels like another one of those situations where Julia is inordinately favoring
interactive usage (REPL, notebooks...) over "batch" usage, even though both of
these usages have their place in practice.

### Surprising performance

Adding to my previous remarks about performance ergonomics, I must say that I
found Julia performance hard to reason about, in the sense that two seemingly
equivalent formulations of a given operation that any competent compiler for a
statically typed and AoT-compiled language would have optimized into the same
machine code would exhibit dramatic run-time performance differences in Julia.

Sometimes, it would just be a matter of dynamic typing not being a zero-cost
abstraction, in the sense that its mere existence slows down your code even when
you're not using it. Many times, however, that was not sufficient to explain
the observed behavior. The compiler had all the information available to
optimize out one particular layer of abstraction, but would for some reason fail
to leverage this knowledge.

One particularly puzzling area was inlining. The Julia compiler turned out to be
surprisingly shy about automatically inlining small functions, while manual
`@inline` hints would have strange side-effects like breaking some heap
allocation elision somewhere, leading to more allocations and a net slowdown.
I'll come back to this in a dedicated section.

### JIT vs profiling tradeoff

Continuing on performance ergonomics, and somewhat complementary to the previous
remark about JIT latency, I found that Julia's highly dynamic operation made it
sometimes hard to distinguish between performance problems in the program under
study and the Julia implementation.

For example, I was persuaded for a while that I was dealing with a JIT corner
case as my program started taking minutes to run, and all stack traces when
manually interrupting it with Ctl+C pointed into something that looked like a
compiler implementation from the file paths. But it later turned out that this
was actually an accidental dynamic typing issue in my code.

Conversely, to get performance profiles that accurately represented the CPU
usage of JITted code, as opposed to that of the JIT implementation, I had to
implement weird JIT warmup code paths that polluted several of my methods with
a `jit_warmup` optional parameter that disabled side-effectful operations
(otherwise stdout printouts and file writes would occur twice).

The fact that the performance contributions of the compiler and the application
code are intermixed and hard to separate like this is arguably a weak point of
the JIT compilation model, from this performance analyst's point of view.

### `@time`

I was pleasantly surprised by `@time` as a very nice quick and dirty performance
analysis tool:

- It's very easy to use as it does not disturb the target expression, which
  continues to operate as before (which cannot be said of, say, `@profile`).
- Its run-time overhead is reasonable.
- It provides useful information about allocations, which are a common
  performance problem in Julia as many language constructs will implicitly lead
  to heap allocations, and GC cycles aren't free.

### Perf investigation tools other than `@profile`

I have the impression that experienced Julia performance engineers are very
quick to study the compiler's various intermediate representations
(`@code_lowered`, `@code_typed`, `@code_warntype`, `@code_llvm` and
`@code_native`). I'm not sure what to think about that:

- On one hand, it's great that this information is so readily available, and it
  can be helpful as a micro-optimization tool or when someone is working on the
  Julia implementation itself.
- On the other hand, the fact that it only concerns a single function call at
  a time makes life harder for those who aren't using the REPL (I know, this is
  becoming a recurring theme), because you need to pay 10s of JITting time on
  every function of a call stack, instead of being able to just dump the
  representation for the entire program and quickly grep through the various
  functions, as is possible with most AoT compilers.
- It also feels like for performance debugging purposes, it could be useful to
  have one layer of intermediate representation between `@code_warntype` (which
  is very easy to use even by complete Julia beginners, but is very close to the
  original code except for type annotations) and `@code_llvm` (which is only
  readable by expert compiler engineers). Something inbetween the moment where
  Julia functions are inlined and the moment where it is translated to LLVM IR.
  This would help debugging inlining issues more quickly, for example.
- And it feels sad that one needs to resort to those very manual tools for now,
  instead of being able to leverage more automated performance analysis tools
  like MAQAO, perf stat (too polluted by JIT right now), etc.
    * I know about `@profile`, but... I'll give it its own part in this report.

Before you ask, I tried out [Traceur](https://github.com/MikeInnes/Traceur.jl),
and found it unpleasant to use in practice because...

- It is very slow compared to direct code execution, and thus does not seem
  applicable beyond very narrow micro-optimization use cases or scenarios where
  you can afford to tune down loop iteration counts.
    * That's sad, because intuitively such tracing use cases sound like a
      situation where JITting could come in handy and allow generating
      reasonably fast instrumented code.
    * Tuning down loop iteration counts can also hide problems which are rare
      but highly impactful, such as occasional generation of denormals.
- It does not handle loops very well, reporting the same errors many times
  instead of deduplicating its output.
- It does not clearly state where the reported errors are coming from, figuring
  it out requires hours of printf instrumentation or trial and error.

Overall, it did not feel like it saved me much time over diving through the
codebase with the performance tips at hand and manually instrumenting the code,
which is bad news for a tool that tries to automate all this manual and
error-prone busywork.

For comparison, Julia's `--track-allocation` mode proved more useful, because
although its effect on run-time performance is similarly dramatic, it does
produce a clear report of where issues were located in the end, in the form of
`.mem` files. So you can just grab a coffee, let it run unattended, and check
out what it found in the end.

### Profiling support

I was deeply disappointed by Julia's profiling support, and am surprised that
this is not a more regular topic of concern for fellow performance engineers who
joined the Julia bandwagon. Let me explain why:

- From what I could observe, `@profile` is the only tool that you can use
  without a custom Julia build. I won't bother rebuilding my Linux distro's
  Julia package just to be able to use my usual tools like perf, and I can tell
  you that the researcher colleagues for the benefit of whom I'm studying Julia
  do not have the time and skills to do so either.
    * Just enable it by default. If you do but Linux distros don't, work with
      them to have their defaults change.
- In its default configuration, `Profile.print()` produces kilometers of
  very poorly organized output by default.
    * It takes a lot of trial-and-error optional parameter tuning to make it
      output something reasonably short and readable by a human, and even then
      the result isn't great.
    * Lines of code are not the right default granularity, there are too many of
      them while function-based granularity is usually enough.
        - Further, in Julia's highly polymorphic design, LoC sorting hides
          important "type of function arguments" information.
        - Yet, at this point in time, function-based profiles are not even
          provided as an option...
    * Lines of code are a bad default sorting keys, people looking at a profile
      are usually most interested in "where do I spend most of my time" so
      sample counts should be the default sorting key at a given stack depth.
        - Yet at this point in time, sample count-based sorting is not even
          provided as an option for hierarchical profiles.
    * Lack of interactivity means that `Profile.print` must throw the entire
      hierarchy at you, complicating high-level analysis which should be the
      first thing that one goes for. `maxdepth` helps, but again setting it up
      is trial and error.
    * Single-space indents are not enough to visually separate stack frames
      well, especially when one tries to figure out which parts of the profile
      occur at a certain stack depth.
    * Absolute sample counts, although useful as an indicator of statistical
      significance, are meaningless per se. Yet there is no option to highlight
      the relative importance of a given part of a profile, either by displaying
      sample percentages or even just color-highlighting sample counts.
    * It does not do a good job at flattening deep call chains where nothing is
      going on (same sample count at each layer of the call stack), as e.g.
      `perf report` does.
    * It does not go a good job at highlighting areas of large self overhead
      (frame at stack depth N has many more samples than the sum of its children
      at stack depth N+1), which are an important target for optimization.
    * It does not go a good job at highlighting which functions were inlined and
      which weren't, even though the julia JIT's highly conservative inlining
      heuristics could have a major impact on performance when a leaf function
      is called very frequently.
- I tried [ProfileView](https://github.com/timholy/ProfileView.jl), and it would
  just hang forever on my machine. But from a look at the README screenshots, it
  does look quite crude and alpha-quality anyhow...

To summarize, I think the measurement back-end of Julia's built-in profiler is
sound, but its user interface needs serious work, and it would be better if
other popular profilers like `perf` worked on Julia programs _by default_.

### Inlining weirdness

From a look at things like `@code_native`, there is surprisingly little inlining
going on in the Julia version. Given how critical inlining is to the performance
of this code (short integration loop iterations, opportunities for
specialization and check elision...), I wouldn't be surprised if it accounted
for a significant fraction of the performance difference wrt the Rust version.

My current hypothesis is that this small degree of inlining might be explained
by the very dynamic nature of Julia's JIT, which may only feed LLVMs with tiny
pieces of code at a time and may not always manage to infer which version of a
given method will be called in the future.

More generally, I could also see this dynamic nature interacting negatively with
the compiler's ability to perform large-scale program analysis. But I'm not a
compiler engineer, and that's mostly speculation.

What I did observe, however, is that manual inlining with `@inline` has
surprising ramifications. For example, I once observed a case where inlining a
method led `@code_native` to report an indirect call to an anonymous Julia
method (something like `julia_#2`, rather than a symbol name which was clearly
correlated to a method name in my code base). And in another situation, inlining
led to a method performing new memory allocations, presumably because the
compiler didn't manage to optimize out an allocation which it managed to
optimize before... for some reason.

Those symptoms might be related to the fact that function calls are a code
specialization boundary in Julia, if specialization occurs after inlining rather
than before. But whatever the explanation, the practical consequence is that
Julia inlines poorly by default, and manual inlining can have unintended
side-effects (beyond the usual ICache footprint growth) that lead to net
performance slowdowns.


## Conclusions

So, how well does Julia deliver on its take at numerical computing world
domination, and in which circumstances (if any), would I recommend using it?

To answer this question, let's quickly go over my main areas of interest again:

- Documentation was bad, but not terrible. It is something which I expect any
  specialized computing engineer to survive, but I would be more hesitant to put
  it in the hands of my less technically-minded colleagues.
- Syntax was a mixed bag, with a lot of good ideas and a fair number of equally
  bad ideas. Overall, Julia succeeds at its somewhat unique goal of looking more
  like math than any other language I've tried, but it does feels carelessly
  designed in other places compared with other contemporary programming
  languages. Nothing feels downright terrible though.
- I feel that programming in the large is a major weak point of Julia as of
  version 1.3, and until this is improved I would strongly hesitate to recommend
  use of this language in any >10kLOC project. Most likely this area hasn't
  received much work because small-scale interactive use cases like notebooks
  and REPL have been prioritized so far.
- The Julia type system was satisfactory overall, but felt weak in my personal
  pet area of automated program correctness checking. I would personally
  recommend future language evolutions to provide more static type checking and
  linting features. They may be optional if not anyone is on board with the
  idea, but to those for whom they matter, they _really_ matter.
- The standard library and wider ecosystem felt solid overall, aside from a few
  bits of sadness here and there. As with the core language, documentation
  quality is a problem that could use more attention.
- Performance ergonomics was the area where Julia disappointed me the most. Good
  performance is often presented as a killer feature of the Julia language,
  giving it an advantage over the more mature Python ecosystem and possibly
  eliminating the need for dual-language applications (e.g. C++ backend +
  Python frontend). However, in practice...
    * Julia's dynamic features lead, even when unused, to a huge amount of
      performance pitfalls that a programmer must constantly bear in mind in
      order to keep their program's performance competitive with that of a
      program written in an AoT-compiled and statically typed language.
    * The Julia JIT is very slow, especially when used non-interactively, which
      makes any file-based development workflow extremely unpleasant. A code
      cache that is retained across program runs would be much appreciated here.
    * Julia performance depends on program structure and design decisions much
      more strongly than in other performance-minded languages, which means that
      it is hard to write Julia code which is simultaneously idiomatic,
      straightforward and performant.
    * The performance analysis tools provided by Julia are, with the exception
      of `@time`, very unpleasant to use and make sense of. Code-based tools
      would benefit from a recursive mode, and profiling tools would benefit
      from a deep user interface rework.

So overall...

- I think that Julia lives up well to its promise of providing much better
  performance than Python given equal development effort, but not to its promise
  of being competitive with "static" languages on the performance front to the
  point of making use of said languages unnecessary.
- I can see why someone would use Julia over Python in interactive use cases
  like notebooks and REPL, but the language manual and ecosystem documentation
  needs to improve significantly before I can personally recommend the language
  to colleagues for this purpose.
- I think that Julia, in its current state at least, is not yet a solid
  contender to statically typed and AoT-compiled languages for writing
  applications, large scale frameworks, and tacking performance-sensitive use
  cases. This is due to its low abilities in areas such as programming in the
  large, performance ergonomics, and automated error-checking.
