module mod_socket
    !socket
    INTEGER:: sock_socket, sock_inet, sock_port        ! socket ID & address of the socket
    CHARACTER(LEN=1024) :: sock_host
    INTEGER, PARAMETER  :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
    CHARACTER(LEN=60)   :: sock_extra_string="                                                            "
    real(8)             :: sock_ecutwf(2)
end module
