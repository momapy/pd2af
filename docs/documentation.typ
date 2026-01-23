#import "@preview/rubber-article:0.5.0": *

#show: article.with(
  cols: none, // Tip: use #colbreak() instead of #pagebreak() to seamlessly toggle columns
  eq-chapterwise: true,
  eq-numbering: "(1.1)",
  header-display: true,
  header-title: "From Process Description to Activity Flow",
  lang: "en",
  page-margins: 1.75in,
  page-paper: "a4",
)

#maketitle(
  title: "From Process Description to Activity Flow",
  authors: ("Adrien Rougny",),
  date: datetime.today().display("[day]/[month]/[year]"),
)

= Introduction

= Method

We first find which entity pools or phenotypes are active (i.e., perform an activity).
We then find influences between activities.

== Finding active entity pools

An entity pool is active _iff_:
- it has a state variable with value "active";
- it has the active state (CellDesigner);
- it is a phenotype;
- it is the source of a modulation;
- it forms a complex with an entity pool via an association; this entity pool is active; and the complex is not active.

== Finding influences between activities

An activity $A$ influences an activity $B$ _iff_:
- the entity pool or phenotype performing $A$ modulates the entity pool or phenotype performing $B$;
- the entity pool performing $A$ modulates one or more processes involved in the production of the entity pool performing $B$;
- the entity pool performing $A$ modulates one or more processes involved in the consumption of the entity pool performing $B$.
