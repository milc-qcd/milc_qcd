/**************** get_qmp_node_number.c ****************************/
/* MIMD version 7 */

/*--------------------------------------------------------------------------*/
/* Get the QMP node number for this machine */
/* C. DeTar 1/30/05 */
/* Within the job script, run this code on all nodes */
/* Prints the qmp logical node number in the form nnnn with leading zeros */
/* WARNING: this can be done only once in any mpirun execution */

#include <stdio.h>
#include <qmp.h>

int
main (int argc, char** argv)
{
  QMP_status_t status;
  int this_node;
  QMP_thread_level_t req, prv;

  /* Start QMP */
  req = QMP_THREAD_SINGLE;
  status = QMP_init_msg_passing (&argc, &argv, req, &prv);

  if (status != QMP_SUCCESS) {
    QMP_error ("QMP_init failed: %s\n", QMP_error_string(status));
    QMP_abort(1);
  }

  /* Get my logical node number */
  this_node = QMP_get_node_number();

  /* Print the result */
  printf("%04d",this_node);

  /* Quit */
  QMP_finalize_msg_passing ();

  return 0;
}
