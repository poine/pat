import rclpy


def main(args=None):
    rclpy.init(args=args)
    wp_publisher = WPPublisher()
    rclpy.spin(wp_publisher)
    wp_publisher.destroy_node()
    rclpy.shutdown()


if __name__ == '__main__':
    main()
