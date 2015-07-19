package utilities;

public interface Flag<T> {
	
	public String getName();
	
	public void setValue(T value);
	
	public T getValue();
	
	public String getDescription();
}
