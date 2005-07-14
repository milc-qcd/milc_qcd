void
restart_gather_from_temp_special(
  void *field,		/* which field? Pointer returned by malloc() */
  int size,		/* size in bytes of the field (eg sizeof(su3_vector))*/
  int index,		/* direction to gather from. eg XUP - index into
			   neighbor tables */
  int parity,		/* parity of sites whose neighbors we gather.
			   one of EVEN, ODD or EVENANDODD. */
  char ** dest,		/* one of the vectors of pointers */
  msg_tag *mtag)          /* previously returned by start_gather */
{
  int i,j;
  site *s;
  msg_t *mbuf;
  gather_t *gt;         /* pointer to current gather */

  //if(mynode()==0) {printf("begin declare: %i\n",index);fflush(stdout);}
  PRINT("restart_gather_from_temp_special\n");
  if(mtag->nsends!=0) mbuf = mtag->send_msgs;
  else mbuf = NULL;

  /* sanity checks for improper usage */
  if(mbuf!=NULL) {
    if(size!=mbuf->gmem->stride) {
      printf("error: wrong stride in restart gather\n");
      terminate(1);
    }
    if(size!=mbuf->gmem->size) {
      printf("error: wrong size in restart gather\n");
      terminate(1);
    }
  }
  gt = &gather_array[index];

  /* set pointers in sites whose neighbors are on this node.  (If all
     neighbors are on this node, this is the only thing done.) */
  if(parity==EVENANDODD) {
    FORALLSITES(j,s){ if(gt->neighbor[j] != NOWHERE){
      dest[j] = field + gt->neighbor[j]*size;
    }}
  } else {
    FORSOMEPARITY(j,s,parity){ if(gt->neighbor[j] != NOWHERE){
      dest[j] = field + gt->neighbor[j]*size;
    }}
  }

  do_gather(mtag);
}

int *neighbor(int dir)
{
  return (gather_array[dir].neighbor);
}

